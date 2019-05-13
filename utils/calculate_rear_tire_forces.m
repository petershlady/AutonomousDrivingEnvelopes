function [Ca_lin, F] = calculate_rear_tire_forces(alpha_prev, P)
% Function to calculate the rear tire forces for each value of k. As the
% rear tire slip angles are not known a priori, it uses the slip angles
% from the previous solution and a linearization of the cornering stiffness
% at these points.
% 
% Inputs:
%   alpha_prev:         k-vector of slip angles at previous time step
%   P:                  parameter struct
% 
% Outputs:
%   Ca_lin:             k-vector containing tire cornering stiffnesses
%                       linearized at each previous slip angle
%   F:                  k-vector containing rear tire force from previous
%                       slip angles
% 
% Usage:
%   [Ca_lin, F] = calculate_rear_tire_forces(alpha_prev, P);
% 
% History:
%   Peter Schleede, 4/17/19 - Initial version
%   Peter Schleede, 4/18/19 - Finished linearization
%   Peter Schleede, 4/21/19 - Fixed an issue with linearization sign
%   Peter Schleede, 5/08/19 - Modified output arguments to allow
%                             integration into dynamics

% get force at previous time steps
F = f_tire(alpha_prev, P.veh.tire_mode, P);

% get linearized cornering stiffnesses (finite differencing)
h = 0.001;
forward = f_tire(alpha_prev + h, P.veh.tire_mode, P);
backward = f_tire(alpha_prev - h, P.veh.tire_mode, P);
Ca_lin = (backward - forward) / (2*h);

end