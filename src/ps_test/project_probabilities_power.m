function [p_out, info] = project_probabilities_power(p_in, x, P_avg, p_min, opts)
%PROJECT_PROBABILITIES_POWER  Project probabilities onto PS constraints.
%
%   p_out = project_probabilities_power(p_in, x, P_avg, p_min)
%
%   Projects an arbitrary vector p_in onto the feasible set:
%       p_i >= p_min
%       sum_i p_i = 1
%       sum_i p_i x_i = P_avg
%
%   The projection solves:
%       minimize ||p - p_in||^2
%       subject to the constraints above.
%
%   The default implementation uses a fast active-set projection for the
%   equality constraints + lower bounds. If it fails, it falls back to
%   fmincon when available.
%
%   Inputs:
%       p_in   - candidate probability vector, Mx1 or 1xM
%       x      - fixed constellation levels, Mx1 or 1xM
%       P_avg  - average optical power constraint
%       p_min  - lower bound for probabilities, e.g. 1e-8
%       opts   - optional struct:
%                  opts.tol        default 1e-10
%                  opts.maxActive  default M
%                  opts.verbose    default false
%
%   Outputs:
%       p_out  - feasible probability vector, Mx1
%       info   - diagnostic struct

    if nargin < 4 || isempty(p_min), p_min = 0; end
    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'tol'),       opts.tol = 1e-10; end
    if ~isfield(opts, 'verbose'),   opts.verbose = false; end

    x = real(x(:));
    M = numel(x);
    p_in = real(p_in(:));

    if numel(p_in) ~= M
        error('project_probabilities_power:SizeMismatch', ...
            'p_in and x must have the same length.');
    end

    if ~all(isfinite(x)) || ~isfinite(P_avg)
        error('project_probabilities_power:BadInput', ...
            'x and P_avg must be finite.');
    end

    p_in(~isfinite(p_in)) = 1/M;
    p_min = max(0, p_min);

    info = struct();
    info.method = 'active-set';
    info.success = false;
    info.message = '';
    info.iterations = 0;

    % Quick feasibility check.
    if M * p_min > 1 + opts.tol
        error('project_probabilities_power:InfeasiblePmin', ...
            'M*p_min must be <= 1. M*p_min = %.6g', M*p_min);
    end

    mass_free = 1 - M*p_min;
    mean_free_target = P_avg - p_min * sum(x);

    if mass_free < opts.tol
        p_fixed = p_min * ones(M,1);
        if abs(sum(p_fixed)-1) < 1e-8 && abs(x.'*p_fixed-P_avg) < 1e-8
            p_out = p_fixed;
            info.success = true;
            info.message = 'All probabilities fixed at p_min.';
            return;
        else
            error('project_probabilities_power:InfeasiblePmin', ...
                'p_min leaves no free probability mass but power constraint is not satisfied.');
        end
    end

    target_avg_free = mean_free_target / mass_free;
    if target_avg_free < min(x) - 1e-8 || target_avg_free > max(x) + 1e-8
        error('project_probabilities_power:InfeasiblePower', ...
            ['No feasible probability vector exists for this x, P_avg and p_min. ', ...
             'Required free average %.8g is outside [%.8g, %.8g].'], ...
             target_avg_free, min(x), max(x));
    end

    % Start from a cleaned version of p_in.
    p_ref = p_in;
    if all(p_ref == 0) || ~isfinite(sum(p_ref))
        p_ref = ones(M,1)/M;
    end

    % Active-set projection onto affine constraints and p >= p_min.
    fixed = false(M,1);
    p = zeros(M,1);
    maxIter = M + 2;

    for iter = 1:maxIter
        info.iterations = iter;
        free = ~fixed;
        F = find(free);
        K = numel(F);

        p(fixed) = p_min;

        b = [1 - sum(p(fixed));
             P_avg - x(fixed).'*p(fixed)];

        if K == 0
            break;
        end

        AF = [ones(1,K); x(F).'];
        pF_ref = p_ref(F);

        % Equality-constrained Euclidean projection:
        % pF = pF_ref - AF' * lambda,
        % where lambda solves AF*AF'*lambda = AF*pF_ref - b.
        G = AF * AF.';
        rhs = AF * pF_ref - b;

        if rcond(G) < 1e-12
            % Degenerate case, for example all remaining x values identical.
            % Use pseudo-inverse and later verify feasibility.
            lambda = pinv(G) * rhs;
        else
            lambda = G \ rhs;
        end

        p(F) = pF_ref - AF.' * lambda;

        viol = find(free & (p < p_min - opts.tol));
        if isempty(viol)
            % Numerical cleanup.
            p(p < p_min) = p_min;

            % One final equality correction on all non-bound entries.
            free2 = p > p_min + 100*eps;
            if any(free2)
                F2 = find(free2);
                A2 = [ones(1,numel(F2)); x(F2).'];
                b2 = [1 - sum(p(~free2)); P_avg - x(~free2).'*p(~free2)];
                G2 = A2*A2.';
                rhs2 = A2*p(F2) - b2;
                if rcond(G2) >= 1e-12
                    p(F2) = p(F2) - A2.'*(G2\rhs2);
                end
            end

            if is_feasible(p, x, P_avg, p_min, 1e-7)
                p_out = normalize_cleanup(p, x, P_avg, p_min);
                info.success = true;
                info.message = 'Projected by active-set method.';
                return;
            end
        end

        % Fix the most negative violating variable at p_min.
        [~, worstLocal] = min(p(viol));
        fixed(viol(worstLocal)) = true;
    end

    % Fallback: fmincon, useful for rare degenerate active-set cases.
    [p_out, info2] = fallback_fmincon_projection(p_ref, x, P_avg, p_min, opts);
    info.method = info2.method;
    info.success = info2.success;
    info.message = info2.message;
    info.iterations = info.iterations + info2.iterations;
end


function tf = is_feasible(p, x, P_avg, p_min, tol)
    tf = all(isfinite(p)) && ...
         all(p >= p_min - tol) && ...
         abs(sum(p) - 1) <= 10*tol && ...
         abs(x(:).'*p(:) - P_avg) <= 10*tol;
end


function p = normalize_cleanup(p, x, P_avg, p_min)
    % Small numerical cleanup. If the active-set result is already feasible,
    % keep it essentially unchanged.
    p = real(p(:));
    p(abs(p) < 1e-14) = 0;
    p(p < p_min) = p_min;

    % Do not do simple p/sum(p), because it may break the power constraint.
    if abs(sum(p)-1) > 1e-9 || abs(x(:).'*p - P_avg) > 1e-9
        try
            [p2, info] = fallback_fmincon_projection(p, x, P_avg, p_min, struct('tol',1e-12,'verbose',false));
            if info.success
                p = p2;
            end
        catch
            % Keep the current value; caller can inspect constraints.
        end
    end
end


function [p_out, info] = fallback_fmincon_projection(p_ref, x, P_avg, p_min, opts)
    info = struct('method','fmincon-fallback', 'success',false, ...
                  'message','', 'iterations',0);

    M = numel(x);
    Aeq = [ones(1,M); x(:).'];
    beq = [1; P_avg];
    lb = p_min * ones(M,1);
    ub = ones(M,1);

    % Build a simple feasible initial point by solving a linear distribution
    % between two constellation points that bracket P_avg, then add p_min.
    p0 = feasible_initial_probability(x, P_avg, p_min);

    if exist('fmincon', 'file') ~= 2
        error('project_probabilities_power:ProjectionFailed', ...
            ['Active-set projection failed and fmincon is unavailable. ', ...
             'Try p_min=0 or check feasibility.']);
    end

    objective = @(p) sum((p(:) - p_ref(:)).^2);
    options = optimoptions('fmincon', ...
        'Display', 'none', ...
        'Algorithm', 'sqp', ...
        'MaxIterations', 100, ...
        'MaxFunctionEvaluations', 2000, ...
        'OptimalityTolerance', 1e-12, ...
        'ConstraintTolerance', 1e-12, ...
        'StepTolerance', 1e-12);

    [p_sol, ~, exitflag, output] = fmincon(objective, p0, [], [], Aeq, beq, lb, ub, [], options);
    info.iterations = output.iterations;
    info.success = exitflag > 0 && is_feasible(p_sol, x, P_avg, p_min, 1e-7);
    info.message = sprintf('fmincon exitflag=%d', exitflag);

    if ~info.success
        error('project_probabilities_power:FminconFailed', ...
            'Probability projection failed. %s', info.message);
    end

    p_out = p_sol(:);
end


function p0 = feasible_initial_probability(x, P_avg, p_min)
    x = x(:);
    M = numel(x);
    p0 = p_min * ones(M,1);

    mass_free = 1 - M*p_min;
    mean_free_target = P_avg - p_min*sum(x);

    if mass_free <= 0
        return;
    end

    target = mean_free_target / mass_free;

    [xs, idx] = sort(x, 'ascend');

    if target <= xs(1)
        p0(idx(1)) = p0(idx(1)) + mass_free;
        return;
    elseif target >= xs(end)
        p0(idx(end)) = p0(idx(end)) + mass_free;
        return;
    end

    hi = find(xs >= target, 1, 'first');
    lo = hi - 1;

    if abs(xs(hi) - xs(lo)) < eps
        p0(idx(lo)) = p0(idx(lo)) + mass_free;
    else
        alpha_hi = (target - xs(lo)) / (xs(hi) - xs(lo));
        alpha_lo = 1 - alpha_hi;
        p0(idx(lo)) = p0(idx(lo)) + mass_free * alpha_lo;
        p0(idx(hi)) = p0(idx(hi)) + mass_free * alpha_hi;
    end
end
