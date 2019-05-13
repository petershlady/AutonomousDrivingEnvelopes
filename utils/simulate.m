function x1 = simulate(x, Fyf, alpha_r, K, dt, f_inv, P)
% Function to simulate one time step and update the state.
% 
% Inputs:
%   x:                  current state
%   Fyf:                front tire force
%   alpha_r:            rear tire slip angle
%   K:                  path curvature
%   dt:                 time step
%   f_inv:              f_tire inverse lookup table
%   P:                  parameter struct
% 
% Ouputs:
%   x1:                 state after time step
% 
% Usage:
%   x1 = simulate(x, Fyf, alpha_r, K, dt, f_inv, P);
% 
% History:
%   Peter Schleede, 4/21/19 - Initial version with linear dynamics
%   Peter Schleede, 4/25/19 - Updated for non-linear dynamics

Fyr = f_tire(alpha_r, 'fiala', P);
alpha_f = calc_f_tire_inv(f_inv, Fyf);
delta = x(1) + P.veh.a*x(2) / P.veh.Ux - alpha_f;
Uy = P.veh.Ux * tan(x(1));

xdot = zeros(5,1);
xdot(1) = (Fyf*cos(delta) + Fyr) / (P.veh.mass*P.veh.Ux) - x(2);
xdot(2) = (P.veh.a*Fyf*cos(delta) - P.veh.b*Fyr) / P.veh.Izz;
xdot(3) = x(2) - P.veh.Ux*K;
xdot(4) = P.veh.Ux*cos(x(3)) - Uy*sin(x(3));
xdot(5) = P.veh.Ux*sin(x(3)) + Uy*cos(x(3));

x1 = x + xdot*dt;

end