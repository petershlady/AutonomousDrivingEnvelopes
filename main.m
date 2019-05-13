% Script to perform the optimization described in "Safe Driving Envelopes
% for Path Tracking in Autonomous Vehicles".
% 
% Author: Peter Schleede
clear

% get parameter struct
run('params');

% initial condition
x0 = [P.x.beta_init, P.x.r_init, P.x.dpsi_init, P.x.s_init, P.x.e_init]';

%% Initialize variables for the problem
% previous alpha_rear will be updated but for now set to zeros
prev_alpha_r    = zeros(P.prob.T_long, 1);

% set things up for stability constraints
% r_min and r_max are based on vehicle properties but beta limits are based
% on r so are introduced as constraints in optimization
alpha_r_peak    = atan2(3*P.veh.mass*9.81*P.veh.mu*P.veh.a, P.veh.Ca*P.veh.L);
r_max           = 9.81 * P.veh.mu / P.veh.Ux;
r_min           = -r_max;

% create f_tire inverse lookup table
f_inv           = create_f_tire_inv(P);

% set fixed curvature for now
% K               = 0*ones(P.prob.T_long+1, 1);
% K               = 0.005*ones(P.prob.T_long+1, 1);
K               = -0.005*ones(P.prob.T_long+1, 1);

%% obstacle
obstacle.s_min  = 100;
obstacle.s_max  = 110;
obstacle.e_min  = -0.5;
obstacle.e_max  = 0.5;

% pre-allocate spaces for data
t_hist          = zeros(1, P.prob.num_steps+1);
x_hist          = zeros(P.prob.num_states, P.prob.num_steps+1);
F_hist          = zeros(1, P.prob.num_steps);
cost_hist       = zeros(5, P.prob.num_steps);
traj_hist       = zeros(P.prob.num_states, P.prob.T_long, P.prob.num_steps);

% set x as x0
x_t = x0;
x_hist(:,1) = x_t;

% variable to help with shaping
Tc = P.prob.T_corr;
Tl = P.prob.T_long;
n_st = P.prob.num_states;
Tcm = Tc - 1;
Tlc = Tl - Tc + 1;

%% Main Loop
%   At each timestep, need to:
%       1. Create time vector
%       2. Create discrete linear system
%       3. Run optimization problem
%       4. Update variables as necessary
tic
for k=1:P.prob.num_steps
    % create time vector
    t = create_time_vector(P);
    if P.prob.t_startup
        P.prob.t_startup = 0;
    end
    
    % get discrete matrices for this time step
    [GsA, GsB1, GsB3, GlA, GlB1, GlB2, GlB3, alpha_r] = ...
        create_discrete_matrices(t, x_t, K, prev_alpha_r, P);
    
    % check for obstacle
    predicted_s = t*P.veh.Ux + x_t(4);
    obstacle_matches = find((predicted_s(Tc+1:end) >= obstacle.s_min) .*...
                         (predicted_s(Tc+1:end) <= obstacle.s_max));
    e_max = P.path.e_max * ones(size(t(Tc+1:end)));
    e_min = P.path.e_min * ones(size(t(Tc+1:end)));
    if ~isempty(obstacle_matches)
        e_min(obstacle_matches) = obstacle.e_max;
    end
    
    %% Optimization problem

    % variables:
    %   x: [beta, r, dpsi, s, e]'
    %       beta:   sideslip
    %       r:      yaw rate
    %       dpsi:   heading deviation
    %       s:      position along path
    %       e:      lateral deviation from path
    %   F_y,f:      forward tire lateral force
    %   v_k = F_yf,k - F_yf,k-1

    % objective: 
    %   sum over k (x_k'*Q*x_k + v_k'*R*v_k + W_veh*sig_veh + W_env*sig_env

    % constraints:
    %   x(k+1)      == A_k*x_k + B_k,1*F_yf,k + B_k,3 for k=1->10
    %   x(k+1)      == A_k*x_k + B_k,1*F_yf,k + B_k,2*F_yf,k+1 + B_k,3 for
    %                       k=11->29
    %   H_veh*x_k   <= G_veh + sig_veh for all k
    %   H_env*x_k   <= G_env + sig_env for k=11->30
    %   F_yf,k      <= F_max for all k
    %   abs(v_k)    <= v_max,k for all k

    % definitions: 
    %   sig_env, sig_veh are non-negative slack variables st the
    %       problem is always feasible
    %   W are weights st the slack variables are always kept small

    % below version of cvx problem is to isolate issues
    cvx_begin quiet
        variables x(5,Tl+1) v(1, Tl) Fyf(1, Tl+1) sig_veh(4) sig_env(2)
        
        % dynamics constraints
        x(:,1) == x_t;
        
        % short time steps
        reshape(x(:,2:Tc), P.prob.num_states*Tcm, 1) == ...
            GsA*reshape(x(:,1:Tcm), P.prob.num_states*(Tcm), 1) ...
            + GsB1*Fyf(1:Tcm)' + GsB3;

        % long time steps
        reshape(x(:,Tc+1:Tl+1), P.prob.num_states*Tlc, 1) == ...
            GlA*reshape(x(:,Tc:Tl), P.prob.num_states*Tlc, 1) ...
            + GlB1*Fyf(Tc:end-1)' + GlB2*Fyf(Tc+1:end)' + GlB3;
        
        % force constraints
        Fyf <= P.con.Fmax;
        Fyf >= P.con.Fmin;
        
        % slew constraints
        v == Fyf(:, 2:end) - Fyf(:, 1:end-1);
        v <= P.con.Vmax;
        v >= P.con.Vmin;
        
        % stability constraints
        x(2, :) <= r_max + sig_veh(1);
        x(2, :) + sig_veh(2) >= r_min;
        x(1, :) <= alpha_r_peak + P.veh.b * x(2, :) / P.veh.Ux  + sig_veh(3);
        x(1, :) + sig_veh(4) >= -alpha_r_peak + P.veh.b * x(2, :) / P.veh.Ux;
        
        % environmental constraints
        x(5, Tc+1:end) <= e_max - (1/2)*P.veh.width - P.path.e_buffer + sig_env(1);
        x(5, Tc+1:end) + sig_env(2) >= e_min + (1/2)*P.veh.width + P.path.e_buffer;
        
        % slack variables
        sig_veh >= 0;
        sig_env >= 0;
        
        % objective
        cost = quad_form(reshape(x(:,2:end), P.prob.num_states*Tl,1), P.opt.Q_tilde)...
               + quad_form(v / P.opt.v_scale, P.opt.R)...
               + P.prob.T_long*P.opt.W_veh' * sig_veh...
               + P.prob.T_long*P.opt.W_env' * sig_env;
        minimize(cost)
    cvx_end
    
    % simulate one time step
    x_t = simulate(x_t, Fyf(1), alpha_r, K(1), P.prob.dt, f_inv, P);
    
    % update previous alpha_r vector
    prev_alpha_r = x(1,2:end) - P.veh.b*x(2,2:end) / P.veh.Ux;
    
    % save observed results
    t_hist(k+1) = k*P.prob.t_short;
    x_hist(:, k+1) = x_t;
    F_hist(:, k) = Fyf(1);
    cost_hist(1, k) = cvx_optval;
    cost_hist(2, k) = sum(x(3,2:end).^2*P.opt.Q_dpsi);
    cost_hist(3, k) = sum(x(5,2:end).^2*P.opt.Q_e);
    cost_hist(4, k) = sum((v/P.opt.v_scale).^2*P.opt.R);
    cost_hist(5, k) = P.prob.T_long*P.opt.W_veh'*sig_veh;
    traj_hist(:, :, k) = x(:, 2:end);
    
    % let user know about progress
    fprintf('Finished step %d/%d, %.3f seconds of simulation time\n',...
                k, P.prob.num_steps, k*P.prob.dt);
    fprintf('CVX status: %s\n', cvx_status);
    
end
toc
%% results visualization
figure(1), clf
hold on
% trajectories
if P.vis.plot_trajs
    for i=1:P.prob.num_steps / P.vis.trajectories
        traj_i = traj_hist(:, :, i*P.vis.trajectories);
        plot(-traj_i(5,:), traj_i(4,:), 'Color', P.vis.colors(3,:))
    end
end
% s and e
plot(zeros(1,P.prob.num_steps+1), x_hist(4,:), 'Color', P.vis.colors(2,:),...
        'LineStyle', '--')
plot(-x_hist(5, :), x_hist(4, :), 'Color', P.vis.colors(1,:), 'LineWidth', 1.5)
% obstacle
plot([-obstacle.e_max, -obstacle.e_min], [obstacle.s_min, obstacle.s_min],...
       'Color', P.vis.colors(4,:), 'LineWidth', 1.5)
plot([-obstacle.e_max, -obstacle.e_min], [obstacle.s_max, obstacle.s_max],...
       'Color', P.vis.colors(4,:), 'LineWidth', 1.5)
plot([-obstacle.e_max, -obstacle.e_max], [obstacle.s_min, obstacle.s_max],...
       'Color', P.vis.colors(4,:), 'LineWidth', 1.5)
plot([-obstacle.e_min, -obstacle.e_min], [obstacle.s_min, obstacle.s_max],...
       'Color', P.vis.colors(4,:), 'LineWidth', 1.5)
xlim([-P.path.e_max,-P.path.e_min])
xlabel('e (m)')
ylabel('s (m)')
title('Error along path')

% heading error
figure(2), clf
plot(t_hist, x_hist(3, :), 'Color', P.vis.colors(1,:))
xlabel('time')
ylabel('\Delta\Psi')
title('Heading error')

% objective
figure(3), clf
plot(t_hist(2:end), cost_hist(1,:), 'Color', P.vis.colors(1,:))
hold on
plot(t_hist(2:end), cost_hist(2,:), 'Color', P.vis.colors(2,:))
plot(t_hist(2:end), cost_hist(3,:), 'Color', P.vis.colors(3,:))
plot(t_hist(2:end), cost_hist(4,:), 'Color', P.vis.colors(4,:))
xlabel('time')
ylabel('Cost')
title('Costs')
legend('Objective cost', 'Heading error cost', 'Lateral error cost', 'Slew rate cost')

% steering angles
delta = rad2deg(x_hist(1, 1:end-1) + P.veh.a*x_hist(2, 1:end-1) /...
            P.veh.Ux - calc_f_tire_inv(f_inv, F_hist));
figure(4), clf
plot(t_hist(1:end-1), delta, 'Color', P.vis.colors(1,:))
xlabel('time')
ylabel('\delta')
title('Steering angles in degrees')

% check stability
beta_max_u = alpha_r_peak + P.veh.b*r_max / P.veh.Ux;
beta_max_l = alpha_r_peak + P.veh.b*r_min / P.veh.Ux;
beta_min_l = -(alpha_r_peak + P.veh.b*r_max / P.veh.Ux);
beta_min_u = -(alpha_r_peak + P.veh.b*r_min / P.veh.Ux);
r = r_min:.01:r_max;
beta = alpha_r_peak + P.veh.b*r ./ P.veh.Ux;
beta_min = -alpha_r_peak +P.veh.b*r ./ P.veh.Ux;

figure(5), clf
hold on
grid on
plot([beta_min_l,beta_max_l], [r_min,r_min], 'b')
plot([beta_min_u,beta_max_u], [r_max,r_max], 'b')
plot(beta, r, 'b')
plot(beta_min, r, 'b')
scatter(x_hist(1, :), x_hist(2, :), 'r', 'x')
xlabel('\beta')
ylabel('r')
title('Stability envelope')