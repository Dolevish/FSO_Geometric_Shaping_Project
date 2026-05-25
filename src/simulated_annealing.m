function out = simulated_annealing(cfg, x0, startIdx, taskLabel)
%SIMULATED_ANNEALING  Single-start SA for 1D constellation.
%
% Progress logging every sa.logEvery iterations with ETA.
% Uses AMI_functions.project_constellation_1D for projection.
% Power constraint enforced inside Enforce_Power_Constraint.m.

if nargin < 4 || isempty(taskLabel)
    taskLabel = sprintf("[SA-%02d]", startIdx);
end

sa = cfg.SA;
tStart = tic;

maxIter      = sa.maxIter;
itersPerTemp = sa.itersPerTemp;

T       = sa.T0;
baseStd = sa.baseStd0;

% Centralized projection
x  = AMI_functions.project_constellation_1D(x0, cfg);
mi = cfg.AMI_Evaluator(x);

bestMI = mi;
x_best = x;

fprintf('%s START | initMI=%.6f | %d iter\n', taskLabel, mi, maxIter);

accCount = 0;
for it = 1:maxIter

    x_prop = x + baseStd * randn(size(x));
    x_prop = AMI_functions.project_constellation_1D(x_prop, cfg);

    mi_prop = cfg.AMI_Evaluator(x_prop);
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

    % --- Progress log ---
    if sa.logEvery > 0 && mod(it, sa.logEvery) == 0
        elapsed = toc(tStart);
        rate    = it / elapsed;
        eta_s   = (maxIter - it) / max(rate, 1);
        fprintf('%s %5d/%d | T=%.2e | MI=%.6f | best=%.6f | %s | ETA %s\n', ...
            taskLabel, it, maxIter, T, mi, bestMI, ...
            fmt_t(elapsed), fmt_t(eta_s));
    end

    % --- Temperature + step adaptation ---
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
    if sec < 60,        s = sprintf('%.0fs', sec);
    elseif sec < 3600,  s = sprintf('%dm%02.0fs', floor(sec/60), mod(sec,60));
    else,               s = sprintf('%dh%02dm', floor(sec/3600), floor(mod(sec,3600)/60));
    end
end
