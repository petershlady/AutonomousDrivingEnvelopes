% Script to store all the parameters used for finding a safe path within
% stability limits and around obstacles. Stores values in a struct that can
% be used elsewhere.
% 
% Author: Peter Schleede

addpath('utils');

%% initial conditions
P.x.beta_init       = 0;
P.x.r_init          = 0;
P.x.dpsi_init       = 0;
P.x.s_init          = 0;
P.x.e_init          = .10;

%% constraints
P.con.Fmax          = 4000;
P.con.Fmin          = -P.con.Fmax;
P.con.Vmax          = 500;
P.con.Vmin          = -P.con.Vmax;

%% problem parameters
P.prob.num_states   = 5;

% time vector used at each step
P.prob.t_short      = 0.05;                 % s
P.prob.t_long       = 0.20;                 % s
P.prob.T_corr       = 10;
P.prob.T_long       = 30;
P.prob.dt           = 0.05;                 % s, actual time step we are using
P.prob.t_startup    = 1;                    % do not ever change in this 
                                            % file, used for time vector 
                                            % initialization
P.prob.num_steps    = 250;                    % number of loop iterations

%% vehicle parameters
P.veh.width         = 1.75;                    % meters
P.veh.mass          = 1200;                 % kg
P.veh.Fz            = 9.81 * P.veh.mass;    % N
P.veh.Ux            = 20;                   % m/s
P.veh.a             = 1.25;                 % m
P.veh.b             = 2.0;                  % m
P.veh.L             = P.veh.a + P.veh.b;    % m
P.veh.mu            = 0.75;
P.veh.Ca            = 80000;
P.veh.Ca_lin        = 110000;
P.veh.tire_mode     = 'fiala';
P.veh.Izz           = (1/12) * P.veh.mass * ((P.veh.a+P.veh.b)^2+2^2);

%% path parameters
P.path.e_buffer     = 0.5;                    % meters
P.path.e_max        = 6;
P.path.e_min        = -2;

%% optimization parameters
% weights for slack variables
P.opt.W_veh         = [1, 1, 1, 1]';
P.opt.W_env         = [3e2, 3e2]';

% quadratic terms
P.opt.Q_dpsi        = 1000;
P.opt.Q_e           = 2;
P.opt.Q             = [0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0;
                       0, 0, P.opt.Q_dpsi, 0, 0;
                       0, 0, 0, 0, 0;
                       0, 0, 0, 0, P.opt.Q_e];
% https://www.mathworks.com/matlabcentral/answers/...
% 324971-forming-a-block-diagonal-matrix-of-one-certain-matrix
Q_tilde             = repmat(P.opt.Q, 1, P.prob.T_long);
Q_tilde             = mat2cell(Q_tilde, size(P.opt.Q,1), repmat(size(P.opt.Q,2),1,P.prob.T_long));
P.opt.Q_tilde       = blkdiag(Q_tilde{:});
P.opt.R             = .1;
P.opt.v_scale       = 100;

%% others
P.vis.colors        = lines(5);
P.vis.plot_trajs    = 1;
P.vis.trajectories  = 1;