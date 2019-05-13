function f_inv = create_f_tire_inv(P)
% Function to set up the lookup table for f_tire inverse. Works with an
% interpolation function for in between.
% 
% Inputs:
%   P:                  parameter struct
% 
% Ouputs:
%   f_inv:              struct containing the lookup table
% 
% Usage:
%   f_inv = create_f_tire_inv(P);
% 
% History:
%   Peter Schleede, 4/25/19 - Initial version

% create alpha vector (intentionally pretty expansive)
alpha = [-0.25:0.001:0.25];

% get tire forces
F = f_tire(alpha, 'fiala', P);

% save results
f_inv.alpha = alpha;
f_inv.F = F;


end