function x_normalized = Enforce_Power_Constraint(x_in, P_target)
% ENFORCE_POWER_CONSTRAINT
% Scales a constellation to strictly meet an Average Power Constraint.
%
% Formula: x_new = x_old * sqrt(P_target / P_current)
%
% Input:
%   x_in     - Input constellation vector (complex/real)
%   P_target - Desired average power (E[|x|^2])
%
% Output:
%   x_normalized - Scaled constellation points

    % 1. Calculate current average power
    % Note: x(:) ensures we treat it as a vector regardless of input shape
    P_current = mean(abs(x_in(:)).^2);
    
    % 2. Safety check for zero-power or empty input
    if P_current < 1e-12
        x_normalized = x_in;
        return;
    end

    % 3. Calculate scaling factor
    scaling_factor = sqrt(P_target / P_current);
    
    % 4. Apply scaling
    x_normalized = x_in * scaling_factor;
end