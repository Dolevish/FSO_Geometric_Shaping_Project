function out = simulated_annealing(cfg, x0, startIdx, taskLabel)
%SIMULATED_ANNEALING  Single-start SA for 1-D constellation optimization.
%
% Commit B adds observability only: the optimization/acceptance logic and
% best-solution definition are unchanged.  The function now records enough
% diagnostics to analyze convergence, repeatability and computational cost.
%
% Commit I.10.1 removes only redundant feasibility assertions immediately
% after project_constellation_1D(). The projector itself already performs the
% same assert before returning, so acceptance logic and feasible candidates
% are unchanged while one duplicate O(M) validation pass per proposal is
% avoided.
%
% Required configuration fields are stored in cfg.SA.  History is sampled
% every cfg.SA.historyEvery iterations (default: itersPerTemp) and always at
% iteration 0 and maxIter.
%
% Returned fields include:
%   x0_projected, initialMI, x_best, bestMI, runtime
%   nAccepted, nImprovingAccepted, nWorseAccepted, nRejected
%   acceptanceRate, startIndex, taskLabel
%   history.iter, history.currentMI, history.bestMI,
%   history.temperature, history.stepStd,
%   history.blockAcceptanceRate, history.cumulativeAcceptanceRate

if nargin < 3 || isempty(startIdx), startIdx = 1; end
if nargin < 4 || isempty(taskLabel)
    taskLabel = sprintf("[SA-%02d]", startIdx);
end
if isstring(taskLabel), taskLabel = char(taskLabel); end

sa = cfg.SA;
tStart = tic;

maxIter      = sa.maxIter;
itersPerTemp = sa.itersPerTemp;

if ~isfield(sa,'historyEvery') || isempty(sa.historyEvery)
    historyEvery = itersPerTemp;
else
    historyEvery = max(1, round(sa.historyEvery));
end

T       = sa.T0;
baseStd = sa.baseStd0;

% Centralized projection: project_constellation_1D validates feasibility
% before it returns, so no duplicate assertion is needed here.
x = AMI_functions.project_constellation_1D(x0, cfg);
mi = cfg.AMI_Evaluator(x);

bestMI = mi;
x_best = x;

x0_projected = x(:);
initialMI    = mi;

fprintf('%s START | initMI=%.6f | %d iter\n', taskLabel, mi, maxIter);

% -------------------------------------------------------------------------
% Diagnostics / observability counters
% -------------------------------------------------------------------------
nAccepted          = 0;
nImprovingAccepted = 0;
nWorseAccepted     = 0;
nRejected          = 0;
blockAccepted      = 0;
lastBlockAccRate   = NaN;

% Preallocate sampled history.  +2 covers iteration 0 and a final sample
% when maxIter is not an exact multiple of historyEvery.
nHistMax = ceil(maxIter / historyEvery) + 2;
hIter    = zeros(nHistMax,1);
hCurrent = nan(nHistMax,1);
hBest    = nan(nHistMax,1);
hT       = nan(nHistMax,1);
hStd     = nan(nHistMax,1);
hBlockAR = nan(nHistMax,1);
hCumAR   = nan(nHistMax,1);
hCount   = 1;

% Initial state (iteration 0).
hIter(hCount)    = 0;
hCurrent(hCount) = mi;
hBest(hCount)    = bestMI;
hT(hCount)       = T;
hStd(hCount)     = baseStd;
hBlockAR(hCount) = NaN;
hCumAR(hCount)   = 0;

for it = 1:maxIter

    x_prop = x + baseStd * randn(size(x));
    % The centralized projector validates feasibility internally.
    x_prop = AMI_functions.project_constellation_1D(x_prop, cfg);

    mi_prop = cfg.AMI_Evaluator(x_prop);
    dE      = mi_prop - mi;

    accepted = false;
    if dE >= 0 || rand() < exp(dE / max(T, eps))
        accepted = true;
        x = x_prop;
        mi = mi_prop;

        nAccepted     = nAccepted + 1;
        blockAccepted = blockAccepted + 1;

        if dE >= 0
            nImprovingAccepted = nImprovingAccepted + 1;
        else
            nWorseAccepted = nWorseAccepted + 1;
        end

        if mi > bestMI
            bestMI = mi;
            x_best = x;
        end
    end

    if ~accepted
        nRejected = nRejected + 1;
    end

    % --- Progress log ---
    if sa.logEvery > 0 && mod(it, sa.logEvery) == 0
        elapsed = toc(tStart);
        rate    = it / max(elapsed, eps);
        eta_s   = (maxIter - it) / max(rate, eps);
        fprintf('%s %5d/%d | T=%.2e | step=%.3e | MI=%.6f | best=%.6f | acc=%.3f | %s | ETA %s\n', ...
            taskLabel, it, maxIter, T, baseStd, mi, bestMI, ...
            nAccepted / it, fmt_t(elapsed), fmt_t(eta_s));
    end

    % --- Temperature + step adaptation ---
    if mod(it, itersPerTemp) == 0
        lastBlockAccRate = blockAccepted / itersPerTemp;
        blockAccepted = 0;

        if lastBlockAccRate > sa.targetAccHi
            baseStd = min(baseStd * sa.baseStdGrow, sa.baseStdMax);
        elseif lastBlockAccRate < sa.targetAccLo
            baseStd = max(baseStd * sa.baseStdShrink, sa.baseStdMin);
        end

        T = max(T * sa.coolingRate, sa.Tf);
    end

    % --- Sample history AFTER the possible temperature/step update ---
    if mod(it, historyEvery) == 0 || it == maxIter
        hCount = hCount + 1;
        hIter(hCount)    = it;
        hCurrent(hCount) = mi;
        hBest(hCount)    = bestMI;
        hT(hCount)       = T;
        hStd(hCount)     = baseStd;
        hBlockAR(hCount) = lastBlockAccRate;
        hCumAR(hCount)   = nAccepted / it;
    end
end

runtime = toc(tStart);

history = struct();
history.iter                     = hIter(1:hCount);
history.currentMI                = hCurrent(1:hCount);
history.bestMI                   = hBest(1:hCount);
history.temperature              = hT(1:hCount);
history.stepStd                  = hStd(1:hCount);
history.blockAcceptanceRate      = hBlockAR(1:hCount);
history.cumulativeAcceptanceRate = hCumAR(1:hCount);
history.sampleEvery              = historyEvery;

out = struct();
out.startIndex          = startIdx;
out.taskLabel           = taskLabel;
out.x0_projected        = x0_projected;
out.initialMI           = initialMI;
out.x_best              = x_best(:);
out.bestMI              = bestMI;
out.finalX               = x(:);
out.finalMI              = mi;
out.runtime              = runtime;
out.nIterations          = maxIter;
out.nAccepted            = nAccepted;
out.nImprovingAccepted   = nImprovingAccepted;
out.nWorseAccepted       = nWorseAccepted;
out.nRejected            = nRejected;
out.acceptanceRate       = nAccepted / max(maxIter,1);
out.history              = history;

fprintf('%s DONE  | best=%.6f | accepted=%d/%d (%.3f) | %s\n', ...
    taskLabel, bestMI, nAccepted, maxIter, out.acceptanceRate, fmt_t(out.runtime));

end

function s = fmt_t(sec)
    if sec < 60,        s = sprintf('%.0fs', sec);
    elseif sec < 3600,  s = sprintf('%dm%02.0fs', floor(sec/60), mod(sec,60));
    else,               s = sprintf('%dh%02dm', floor(sec/3600), floor(mod(sec,3600)/60));
    end
end
