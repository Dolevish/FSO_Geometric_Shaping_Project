function x = project_constellation_weighted_power(x_in, p_fixed, cfg)
%PROJECT_CONSTELLATION_WEIGHTED_POWER
% Projects a 1D IM/DD constellation onto the weighted-power feasible set:
%
%   x_i >= 0
%   x sorted ascending
%   x_{i+1} - x_i >= minGap
%   sum_i p_i * x_i = P_avg
%
% This version preserves the association between probabilities and ordered
% constellation ranks. It uses a monotone residual representation:
%
%   x = b + r,  b_i = (i-1)*minGap,  r_1 <= r_2 <= ... <= r_M, r_i >= 0
%
% Scaling r by a nonnegative scalar preserves both ordering and minGap.

    sa = cfg.SA;
    x  = real(x_in(:));
    p  = p_fixed(:);
    M  = numel(x);

    if numel(p) ~= M
        error('project_constellation_weighted_power:SizeMismatch', ...
            'x and p must have the same length.');
    end

    if ~isfield(sa,'minGap') || isempty(sa.minGap)
        minGap = 0;
    else
        minGap = max(0, sa.minGap);
    end

    P_avg = cfg.P_avg;

    % Probability safety. The PS stage should already return a valid p, but
    % this keeps the projection robust if called independently.
    p(~isfinite(p)) = 0;
    p = max(p, 0);
    if sum(p) <= 0
        p = ones(M,1) / M;
    else
        p = p / sum(p);
    end

    % Baseline that exactly enforces the minimum spacing when residual r=0.
    b = (0:M-1).' * minGap;
    baselinePower = dot(p, b);
    residualPowerTarget = P_avg - baselinePower;

    if residualPowerTarget < -1e-12
        error('project_constellation_weighted_power:InfeasibleMinGap', ...
            ['Infeasible minGap for current p: dot(p,b)=%.12g exceeds ', ...
             'P_avg=%.12g. Reduce minGap or change p_min.'], ...
            baselinePower, P_avg);
    end
    residualPowerTarget = max(0, residualPowerTarget);

    % Step 1: project the raw proposal into a nonnegative sorted direction.
    % We do not permute p; p remains attached to ordered constellation ranks.
    x = sort(max(x, 0), 'ascend');

    % Step 2: express x relative to the minGap baseline.
    residual = x - b;

    % Step 3: enforce residual >= 0 and monotonicity. Then x=b+residual
    % automatically satisfies x_{i+1}-x_i >= minGap.
    residual(1) = max(residual(1), 0);
    for i = 2:M
        residual(i) = max(residual(i), residual(i-1));
    end

    % Step 4: scale only the residual to hit the weighted power exactly.
    residualPower = dot(p, residual);
    if residualPower > 0 && isfinite(residualPower)
        x = b + residual * (residualPowerTarget / residualPower);
    else
        % Deterministic feasible fallback when the proposal collapsed.
        template = (0:M-1).';
        templatePower = dot(p, template);
        if templatePower > 0 && isfinite(templatePower)
            x = b + template * (residualPowerTarget / templatePower);
        else
            % This can only occur for a degenerate p concentrated on the first
            % symbol and residualPowerTarget>0. Use all mass at the last rank.
            template = zeros(M,1);
            template(end) = 1;
            templatePower = dot(p, template);
            if templatePower > 0 && isfinite(templatePower)
                x = b + template * (residualPowerTarget / templatePower);
            else
                x = b;
            end
        end
    end

    % Final numerical cleanup only. Do not sort here; sorting could change the
    % p-weighted power if numerical noise ever created a crossing.
    x = real(x(:));
    x(abs(x) < 1e-14) = 0;
end
