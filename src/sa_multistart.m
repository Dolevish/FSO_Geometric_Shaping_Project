function [out, results] = sa_multistart(cfg, x_baseline, taskLabel)
%SA_MULTISTART  Multi-start SA with optional PARFOR parallelism.
%
% Commit B preserves the existing winner-selection rule (largest FAST AMI)
% and adds observability/reproducibility metadata for every start.
%
% - Forces PROCESS-based parallel pool only.
% - Uses a single master seed (random by default) + substreams per start.
% - Stores initial point, final/best solution, acceptance statistics and
%   sampled convergence history for every restart.
% - Prints the master seed for reproducibility/debug.
%
% IMPORTANT: validated-AMi winner selection is intentionally NOT introduced
% here; that change belongs to Commit C.

if nargin < 3 || isempty(taskLabel)
    taskLabel = "[SA]";
end
if isstring(taskLabel), taskLabel = char(taskLabel); end

sa = cfg.SA;

% ---- Defaults (safe) ----
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

% ---- Choose master seed (client-side only) ----
if isempty(sa.seedInit)
    rng("shuffle");
    masterSeed = randi([0, 2^32-1], 1, 1, "uint32");
else
    masterSeed = uint32(sa.seedInit);
end

fprintf('RNG masterSeed (Threefry) = %u\n', masterSeed);

% ---- Build initial points X0 using the same masterSeed ----
% Separate stream: SA random consumption cannot change restart initialization.
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
    'bestMI', -Inf, ...
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

% ---- Run starts ----
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
    % Serial: same substream scheme as parallel mode.
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

% ---- Pick best (UNCHANGED in Commit B: FAST evaluator only) ----
allMI = [results.bestMI];
[bestMI, idx] = max(allMI);

out = struct();
out.bestStart  = idx;
out.bestMI     = bestMI;
out.bestX      = results(idx).x_best(:);
out.masterSeed = masterSeed;
out.nStarts    = nStarts;
out.totalRuntime = sum([results.runtime]);
out.meanAcceptanceRate = mean([results.acceptanceRate], 'omitnan');
out.bestRun    = results(idx);

% ---- Optional pool cleanup ----
if doPar && poolStartedHere && isfield(sa,'closePoolWhenDone') && sa.closePoolWhenDone
    delete(gcp("nocreate"));
end

end
