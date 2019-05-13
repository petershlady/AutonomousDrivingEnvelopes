function alpha = calc_f_tire_inv(f_inv, F)
% Function to interpolate for f_tire inverse. Just uses interp1 but hides
% the implementation.
% 
% Inputs:
%   f_inv:              struct containing lookup table
%   F:                  force to look up
% 
% Ouputs:
%   alpha:              slip angle corresponding to force
% 
% Usage:
%   alpha = calc_f_tire_inv(f_inv, F);
% 
% History:
%   Peter Schleede, 4/25/19 - Initial version
%   Peter Schleede, 5/03/19 - Updated check for bounds on F for vectors

if any(F < min(f_inv.F)) || any(F > max(f_inv.F))
    ME = MException('Envelopes:calc_f_tire_inv', ...
               'requested force is outside limits of lookup table');
    throw(ME);
end

alpha = interp1(f_inv.F, f_inv.alpha, F);

end