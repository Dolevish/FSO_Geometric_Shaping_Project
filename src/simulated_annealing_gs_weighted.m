function out = simulated_annealing_gs_weighted(cfg, x0, p_fixed, startIdx, taskLabel)
%SIMULATED_ANNEALING_GS_WEIGHTED
% Single-start simulated annealing for GS with fixed non-uniform probabilities.
%
% Difference from simulated_annealing.m:
%   - x is optimized, p is fixed.
%   - Projection enforces the weighted power constraint sum_i p_i*x_i=P_avg.
%
% This is the GS step needed for alternating GS/PS.

if nargin < 5 || isempty(taskLabel)
    taskLabel = sprintf('[WGS-%02d]', startIdx);
end
if isstring(taskLabel), taskLabel = char(taskLabel); end

sa = cfg.SA;
tStart = tic;

maxIter      = sa.maxIter;
itersPerTemp = sa.itersPerTemp;

T       = sa.T0;
baseStd = sa.baseStd0;

p_fixed = p_fixed(:);
x  = project_constellation_weighted_power(x0, p_fixed, cfg);
mi = cfg.GS_Evaluator(x, p_fixed);

bestMI = mi;
x_best = x;

fprintf('%s START | initMI=%.6f | %d iter\n', taskLabel, mi, maxIter);

accCount = 0;
for it = 1:maxIter

    x_prop = x + baseStd * randn(size(x));
    x_prop = project_constellation_weighted_power(x_prop, p_fixed, cfg);

    mi_prop = cfg.GS_Evaluator(x_prop, p_fixed);
    dE      = mi_prop - mi;

    if dE >= 0 || rand() < exp(dE / max(T, eps))
        x = x_prop;
        mi = mi_prop;
        accCount = accCount + 1;

        if mi > bestMI
            bestMI = mi;
            x_best = x;
        end
    end

    if sa.logEvery > 0 && mod(it, sa.logEvery) == 0
        elapsed = toc(tStart);
        rate    = it / max(elapsed, eps);
        eta_s   = (maxIter - it) / max(rate, eps);
        fprintf('%s %5d/%d | T=%.2e | MI=%.6f | best=%.6f | %s | ETA %s\n', ...
            taskLabel, it, maxIter, T, mi, bestMI, fmt_t(elapsed), fmt_t(eta_s));
    end

    if mod(it, itersPerTemp) == 0
        accRate  = accCount / itersPerTemp;
        accCount = 0;

        if accRate > sa.targetAccHi
            baseStd = min(baseStd * sa.baseStdGrow, sa.baseStdMax);
        elseif accRate < sa.targetAccLo
            baseStd = max(baseStd * sa.baseStdShrink, sa.baseStdMin);
        end

        T = max(T * sa.coolingRate, sa.Tf);
    end
end

out = struct();
out.x_best  = x_best(:);
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
