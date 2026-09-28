function [out, results] = sa_multistart(cfg, x_baseline, taskLabel)
%SA_MULTISTART  Multi-start SA with independent post-optimization validation.
%
% Commit C changes the final-selection rule:
%   1) every restart is optimized with the FAST evaluator;
%   2) the best constellation from EVERY restart is independently validated;
%   3) the final winner is selected by VALIDATED AMI, not fast AMI.
%
% The function also preserves the Commit-B observability/reproducibility
% metadata for every start and records both FAST and VALIDATED rankings.
%
% Returned out fields include both explicit rankings:
%   bestStartFast, bestMIFast, bestXFast
%   bestStartValidated, bestMIValidated, bestXValidated
%
% For backward compatibility with callers that use out.bestStart/out.bestMI/
% out.bestX, these aliases now refer to the VALIDATED winner.

if nargin < 3 || isempty(taskLabel)
    taskLabel = "[SA]";
end
if isstring(taskLabel), taskLabel = char(taskLabel); end

sa = cfg.SA;

if ~isfield(cfg,'AMI_Validator') || isempty(cfg.AMI_Validator)
    error('sa_multistart:MissingValidator', ...
        ['cfg.AMI_Validator is required for validated multi-start selection. ' ...
         'Build cfg with build_fso_config().']);
end

% ---- Defaults ----
if ~isfield(sa,'nStarts') || isempty(sa.nStarts), sa.nStarts = 1; end
if ~isfield(sa,'useParallel') || isempty(sa.useParallel), sa.useParallel = false; end
if ~isfield(sa,'seedInit'), sa.seedInit = []; end
if ~isfield(sa,'numWorkers'), sa.numWorkers = []; end

nStarts = sa.nStarts;

fprintf('Running SA multistart: nStarts=%d (parallel=%d)\n', nStarts, sa.useParallel);

% ---- Decide parallel ----
doPar = logical(sa.useParallel);
if doPar
    try
        doPar = license('test','Distrib_Computing_Toolbox') ~= 0;
    catch
        doPar = false;
    end
end

% ---- Start or reuse pool (FORCE Processes) ----
poolStartedHere = false;
p = [];
if doPar
    try
        p = gcp("nocreate");

        if ~isempty(p)
            profNow = string(p.Cluster.Profile);
            if ~strcmpi(profNow, "Processes")
                delete(p);
                p = [];
            end
        end

        if isempty(p)
            if isempty(sa.numWorkers)
                p = parpool("Processes");
            else
                p = parpool("Processes", sa.numWorkers);
            end
            poolStartedHere = true;
        end

        fprintf('Parallel pool active: profile=%s | workers=%d\n', ...
            string(p.Cluster.Profile), p.NumWorkers);
    catch ME
        fprintf('Could not start parallel pool (%s). Falling back to serial.\n', ME.message);
        doPar = false;
        p = [];
    end
end

% ---- Choose master seed ----
if isempty(sa.seedInit)
    rng("shuffle");
    masterSeed = randi([0, 2^32-1], 1, 1, "uint32");
else
    masterSeed = uint32(sa.seedInit);
end

fprintf('RNG masterSeed (Threefry) = %u\n', masterSeed);

% ---- Build restart initial points reproducibly ----
sX0 = RandStream('Threefry','Seed',double(masterSeed));
X0 = cell(nStarts,1);
X0{1} = x_baseline(:);
for k = 2:nStarts
    sX0.Substream = k;
    X0{k} = x_baseline(:) + 0.2*randn(sX0, numel(x_baseline), 1);
end

% ---- Preallocate results ----
emptyHistory = struct( ...
    'iter', [], ...
    'currentMI', [], ...
    'bestMI', [], ...
    'temperature', [], ...
    'stepStd', [], ...
    'blockAcceptanceRate', [], ...
    'cumulativeAcceptanceRate', [], ...
    'sampleEvery', []);

results = repmat(struct( ...
    'startIndex', NaN, ...
    'masterSeed', masterSeed, ...
    'substream', NaN, ...
    'x0_raw', [], ...
    'x0_projected', [], ...
    'initialMI', NaN, ...
    'bestMI', -Inf, ...              % legacy alias: FAST AMI
    'bestMIFast', -Inf, ...
    'bestMIValidated', NaN, ...
    'validationGap', NaN, ...
    'validationRuntime', NaN, ...
    'fastRank', NaN, ...
    'validatedRank', NaN, ...
    'x_best', [], ...
    'finalMI', NaN, ...
    'finalX', [], ...
    'runtime', NaN, ...
    'nIterations', NaN, ...
    'nAccepted', NaN, ...
    'nImprovingAccepted', NaN, ...
    'nWorseAccepted', NaN, ...
    'nRejected', NaN, ...
    'acceptanceRate', NaN, ...
    'history', emptyHistory), nStarts, 1);

% ---- Constant stream for SA randomness ----
if doPar
    constantStream = parallel.pool.Constant(@() ...
        RandStream('Threefry','Seed',double(masterSeed)));
else
    constantStream = [];
end

% ---- Run SA starts ----
if doPar
    parfor k = 1:nStarts
        label = sprintf("%s-%02d", taskLabel, k);

        s = constantStream.Value;
        s.Substream = k;
        RandStream.setGlobalStream(s);

        outk = simulated_annealing(cfg, X0{k}, k, label);

        results(k).startIndex          = k;
        results(k).masterSeed          = masterSeed;
        results(k).substream           = k;
        results(k).x0_raw              = X0{k}(:);
        results(k).x0_projected        = outk.x0_projected(:);
        results(k).initialMI           = outk.initialMI;
        results(k).bestMI              = outk.bestMI;
        results(k).bestMIFast          = outk.bestMI;
        results(k).x_best              = outk.x_best(:);
        results(k).finalMI             = outk.finalMI;
        results(k).finalX              = outk.finalX(:);
        results(k).runtime             = outk.runtime;
        results(k).nIterations         = outk.nIterations;
        results(k).nAccepted           = outk.nAccepted;
        results(k).nImprovingAccepted  = outk.nImprovingAccepted;
        results(k).nWorseAccepted      = outk.nWorseAccepted;
        results(k).nRejected           = outk.nRejected;
        results(k).acceptanceRate      = outk.acceptanceRate;
        results(k).history             = outk.history;
    end
else
    s = RandStream('Threefry','Seed',double(masterSeed));
    RandStream.setGlobalStream(s);

    for k = 1:nStarts
        label = sprintf("%s-%02d", taskLabel, k);

        s.Substream = k;
        RandStream.setGlobalStream(s);

        outk = simulated_annealing(cfg, X0{k}, k, label);

        results(k).startIndex          = k;
        results(k).masterSeed          = masterSeed;
        results(k).substream           = k;
        results(k).x0_raw              = X0{k}(:);
        results(k).x0_projected        = outk.x0_projected(:);
        results(k).initialMI           = outk.initialMI;
        results(k).bestMI              = outk.bestMI;
        results(k).bestMIFast          = outk.bestMI;
        results(k).x_best              = outk.x_best(:);
        results(k).finalMI             = outk.finalMI;
        results(k).finalX              = outk.finalX(:);
        results(k).runtime             = outk.runtime;
        results(k).nIterations         = outk.nIterations;
        results(k).nAccepted           = outk.nAccepted;
        results(k).nImprovingAccepted  = outk.nImprovingAccepted;
        results(k).nWorseAccepted      = outk.nWorseAccepted;
        results(k).nRejected           = outk.nRejected;
        results(k).acceptanceRate      = outk.acceptanceRate;
        results(k).history             = outk.history;
    end
end

% -------------------------------------------------------------------------
% Commit C: independently validate EVERY restart winner.
% Validation is deliberately performed after all stochastic SA runs so the
% validation stage cannot influence RNG streams or the SA trajectories.
% -------------------------------------------------------------------------
fprintf('%s validating %d restart winners with independent evaluator...\n', ...
    taskLabel, nStarts);

validationTotalTimer = tic;
for k = 1:nStarts
    AMI_functions.assert_constellation_feasible(results(k).x_best, cfg, 1e-10);

    tVal = tic;
    miVal = cfg.AMI_Validator(results(k).x_best);
    results(k).validationRuntime = toc(tVal);

    if ~(isscalar(miVal) && isfinite(miVal))
        error('sa_multistart:InvalidValidatedAMI', ...
            'Validation returned an invalid AMI for restart %d.', k);
    end

    results(k).bestMIValidated = miVal;
    results(k).validationGap   = results(k).bestMIFast - miVal;

    fprintf('%s validate %02d/%02d | fast=%.6f | val=%.6f | delta=%+.3e | %.2fs\n', ...
        taskLabel, k, nStarts, results(k).bestMIFast, miVal, ...
        results(k).validationGap, results(k).validationRuntime);
end
validationTotalRuntime = toc(validationTotalTimer);

% ---- Rank restarts under both evaluators ----
allFast = [results.bestMIFast];
allVal  = [results.bestMIValidated];

[~, fastOrder] = sort(allFast, 'descend');
[~, valOrder]  = sort(allVal,  'descend');

for r = 1:nStarts
    results(fastOrder(r)).fastRank     = r;
    results(valOrder(r)).validatedRank = r;
end

[bestMIFast, idxFast] = max(allFast);
[bestMIVal,  idxVal]  = max(allVal);

selectionChanged = idxFast ~= idxVal;
validatedAdvantage = bestMIVal - results(idxFast).bestMIValidated;

absGaps = abs([results.validationGap]);

% ---- Final output ----
out = struct();
out.masterSeed = masterSeed;
out.nStarts    = nStarts;

% Explicit fast-evaluator winner.
out.bestStartFast = idxFast;
out.bestMIFast    = bestMIFast;
out.bestXFast     = results(idxFast).x_best(:);
out.fastWinnerValidatedMI = results(idxFast).bestMIValidated;

% Explicit validated winner: THIS is the final selected solution.
out.bestStartValidated = idxVal;
out.bestMIValidated    = bestMIVal;
out.bestXValidated     = results(idxVal).x_best(:);
out.validatedWinnerFastMI = results(idxVal).bestMIFast;

% Diagnostics quantifying whether fast ranking changed the selected restart.
out.selectionChanged = selectionChanged;
out.validatedAdvantageOverFastWinner = validatedAdvantage;
out.meanAbsFastValidationGap = mean(absGaps, 'omitnan');
out.maxAbsFastValidationGap  = max(absGaps);
out.validationTotalRuntime   = validationTotalRuntime;

% SA runtime diagnostics from Commit B.
out.totalRuntime = sum([results.runtime]) + validationTotalRuntime;
out.totalSARuntime = sum([results.runtime]);
out.meanAcceptanceRate = mean([results.acceptanceRate], 'omitnan');

% Backward-compatible aliases now deliberately point to VALIDATED winner.
out.bestStart = out.bestStartValidated;
out.bestMI    = out.bestMIValidated;
out.bestX     = out.bestXValidated;
out.bestRun   = results(idxVal);

fprintf('%s selection summary | fast winner=%d | validated winner=%d | changed=%d\n', ...
    taskLabel, idxFast, idxVal, selectionChanged);
fprintf('%s final validated AMI=%.6f | fast-winner validated AMI=%.6f | improvement=%.3e\n', ...
    taskLabel, bestMIVal, results(idxFast).bestMIValidated, validatedAdvantage);

% ---- Optional pool cleanup ----
if doPar && poolStartedHere && isfield(sa,'closePoolWhenDone') && sa.closePoolWhenDone
    delete(gcp("nocreate"));
end

end
