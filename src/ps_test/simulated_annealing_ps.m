function out = simulated_annealing_ps(cfg, x_fixed, p0, startIdx, taskLabel)
%SIMULATED_ANNEALING_PS  Single-start simulated annealing for PS.
%
%   Optimizes symbol probabilities p for a fixed constellation x_fixed.
%
%   Constraints enforced after every proposal:
%       p_i >= cfg.PS.p_min
%       sum_i p_i = 1
%       sum_i p_i x_i = cfg.P_avg
%
%   Required fields:
%       cfg.PS
%       cfg.PS_Evaluator = @(x_fixed,p) AMI in bits/symbol
%
%   Output:
%       out.p_best  - best probability vector found
%       out.bestMI  - best fast AMI value found
%       out.runtime - elapsed time

    if nargin < 4 || isempty(startIdx), startIdx = 1; end
    if nargin < 5 || isempty(taskLabel)
        taskLabel = sprintf('[PS-%02d]', startIdx);
    end
    if isstring(taskLabel), taskLabel = char(taskLabel); end

    ps = cfg.PS;
    tStart = tic;

    maxIter      = ps.maxIter;
    itersPerTemp = ps.itersPerTemp;

    T       = ps.T0;
    baseStd = ps.baseStd0;

    x_fixed = x_fixed(:);
    p = project_probabilities_power(p0, x_fixed, cfg.P_avg, ps.p_min);
    mi = cfg.PS_Evaluator(x_fixed, p);

    bestMI = mi;
    p_best = p;

    if ~isfield(ps,'logEvery') || isempty(ps.logEvery)
        ps.logEvery = 500;
    end

    fprintf('%s START | initMI=%.6f | %d iter\n', taskLabel, mi, maxIter);

    accCount = 0;
    for it = 1:maxIter

        p_prop = p + baseStd * randn(size(p));
        p_prop = project_probabilities_power(p_prop, x_fixed, cfg.P_avg, ps.p_min);

        mi_prop = cfg.PS_Evaluator(x_fixed, p_prop);
        dE      = mi_prop - mi;

        if dE >= 0 || rand() < exp(dE / max(T, eps))
            p = p_prop;
            mi = mi_prop;
            accCount = accCount + 1;

            if mi > bestMI
                bestMI = mi;
                p_best = p;
            end
        end

        if ps.logEvery > 0 && mod(it, ps.logEvery) == 0
            elapsed = toc(tStart);
            rate    = it / max(elapsed, eps);
            eta_s   = (maxIter - it) / max(rate, eps);
            fprintf('%s %5d/%d | T=%.2e | MI=%.6f | best=%.6f | %s | ETA %s | min(p)=%.2e\n', ...
                taskLabel, it, maxIter, T, mi, bestMI, ...
                fmt_t(elapsed), fmt_t(eta_s), min(p));
        end

        if mod(it, itersPerTemp) == 0
            accRate  = accCount / itersPerTemp;
            accCount = 0;

            if accRate > ps.targetAccHi
                baseStd = min(baseStd * ps.baseStdGrow, ps.baseStdMax);
            elseif accRate < ps.targetAccLo
                baseStd = max(baseStd * ps.baseStdShrink, ps.baseStdMin);
            end

            T = max(T * ps.coolingRate, ps.Tf);
        end
    end

    out = struct();
    out.p_best  = p_best(:);
    out.bestMI  = bestMI;
    out.runtime = toc(tStart);

    fprintf('%s DONE  | best=%.6f | %s\n', taskLabel, bestMI, fmt_t(out.runtime));
end


function s = fmt_t(sec)
    if sec < 60
        s = sprintf('%.0fs', sec);
    elseif sec < 3600
        s = sprintf('%dm%02.0fs', floor(sec/60), mod(sec,60));
    else
        s = sprintf('%dh%02dm', floor(sec/3600), floor(mod(sec,3600)/60));
    end
end
