function F = f_tire(alpha, mode, P)
% Function to compute the lateral tire force. Uses brush tire model in
% Pacejka (2002).
% 
% Inputs:
%   alpha:      tire slip angle
%   mode:       'linear' or 'fiala'
%   P:          parameter struct
% 
% Outputs:
%   F:          lateral tire force
% 
% Usage:
%   F = f_tire(alpha, P);
% 
% History:
%   Peter Schleede, 4/2/19 - Initial version
%   Peter Schleede, 4/17/19 - Added vector support
%   Peter Schleede, 4/18/19 - Added linear option

Ca = P.veh.Ca;
mu = P.veh.mu;
Fz = P.veh.Fz;

if strcmp(mode, 'linear')
    F = -Ca * alpha;
elseif strcmp(mode, 'fiala')
    alpha_sl = atan2(3*mu*Fz, Ca);
    k = length(alpha);
    F = zeros(size(alpha));

    for i=1:k
        ta = tan(alpha(i));
        if abs(alpha(i)) < alpha_sl
            F(i) = -Ca*ta + ...
                    Ca^2 / (3*mu*Fz)*abs(ta)*ta - ...
                    Ca^3 / (27*mu^2*Fz^2)*ta^3;
        else
            F(i) = -mu * Fz * sign(alpha(i));
        end
    end
else
    ME = MException('Envelopes:f_tire', ...
               'tire force mode must be either ''linear'' or ''fiala''');
    throw(ME);
end