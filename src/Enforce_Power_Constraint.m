function x = Enforce_Power_Constraint(x, P_avg, mode)
%ENFORCE_POWER_CONSTRAINT
% mode:
%  "meansquare" -> mean(x.^2)=P_avg (baseband electrical power)
%  "mean"       -> mean(x)=P_avg    (IM/DD average intensity)

if nargin < 3, mode = "meansquare"; end
x = x(:);

switch string(mode)
    case "meansquare"
        p = mean(x.^2);
        if ~(isfinite(p) && p > 0)
            error("Enforce_Power_Constraint:BadPower","mean(x.^2) invalid: %g", p);
        end
        x = x * sqrt(P_avg / p);

    case "mean"
        m = mean(x);
        if ~(isfinite(m) && m > 0)
            error("Enforce_Power_Constraint:BadPower","mean(x) invalid: %g", m);
        end
        x = x * (P_avg / m);

    otherwise
        error("Enforce_Power_Constraint:BadMode","Unknown mode: %s", mode);
end
end
