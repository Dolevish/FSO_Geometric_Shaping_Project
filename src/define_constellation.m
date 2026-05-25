function x = define_constellation(M, P_avg, imdd_mode, powerConstraint)
if nargin < 3, imdd_mode = false; end
if nargin < 4, powerConstraint = "meansquare"; end

if imdd_mode
    % IM/DD intensity levels (nonnegative)
    x = 0:(M-1);
else
    % Bipolar PAM
    x = -(M-1):2:(M-1);
end

x = Enforce_Power_Constraint(x(:), P_avg, powerConstraint);
x = sort(x(:)).';
end
