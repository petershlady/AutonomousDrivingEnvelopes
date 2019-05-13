function t_vect = create_time_vector(P)
% Function to create the time vector for use in MPC model. Follows method
% described in Erlien (2013). At initialization, t_corr can just be
% t_short. However, from then on, it adjusts to make sure that the long
% time steps are the same between control loops
% 
% Inputs:
%   P:          parameter struct
% 
% Outputs:
%   t_vect:     time vector, with t_corr appropriately chosen
% 
% Usage:
%   t_vect = create_time_vector(P);
% 
% History:
%   4/2/19 - Initial version
%   4/19/19 - Modified to have dt as a parameter
%   4/24/19 - Modified to allow for fixed time step
%   5/12/19 - Removed fixed time step

persistent t_corr_location;

% vector of short times is always the same
short_t_vect = P.prob.t_short:P.prob.t_short:(P.prob.T_corr-1)*P.prob.t_short;

% handle moving t_corr around
if P.prob.t_startup
    % figure out where t_corr is for this and where the fixed time point
    % should be
    t_corr = P.prob.T_corr * P.prob.t_short;

    t_corr_location = t_corr + P.prob.t_long;

    % create the full vector
    middle_t_vect = [short_t_vect, t_corr];
    ending_time = (P.prob.T_long-P.prob.T_corr)*P.prob.t_long;
    t_vect = [middle_t_vect, (P.prob.t_long:P.prob.t_long:ending_time)+...
                    middle_t_vect(end)];
else
    % update t_corr_location
    t_corr_location = t_corr_location - P.prob.dt;
    if (t_corr_location-short_t_vect(end)) < P.prob.t_short
        t_corr_location = t_corr_location + P.prob.t_long;
    end

    t_corr = t_corr_location;

    % create the full vector
    middle_t_vect = [short_t_vect, t_corr];
    ending_time = (P.prob.T_long-P.prob.T_corr)*P.prob.t_long;
    t_vect = [middle_t_vect, (P.prob.t_long:P.prob.t_long:ending_time)+...
                    middle_t_vect(end)];
end

t_vect = [0, t_vect];

end