function [out, results] = gs_weighted_multistart(cfg, x_baseline, p_fixed, taskLabel)
%GS_WEIGHTED_MULTISTART
% Multi-start SA for weighted-power GS with fixed probabilities p_fixed.
%
% This is the alternating-GS counterpart of sa_multistart.m.

if nargin < 4 || isempty(taskLabel)
    taskLabel = '[WGS]';
end
if isstring(taskLabel), taskLabel = char(taskLabel); end

sa = cfg.SA;
if ~isfield(sa,'nStarts') || isempty(sa.nStarts), sa.nStarts = 1; end
if ~isfield(sa,'useParallel') || isempty(sa.useParallel), sa.useParallel = false; end
if ~isfield(sa,'seedInit'), sa.seedInit = []; end
if ~isfield(sa,'numWorkers'), sa.numWorkers = []; end

nStarts = sa.nStarts;
fprintf('Running weighted-GS multistart: nStarts=%d (parallel=%d)\n', nStarts, sa.useParallel);

% Decide parallel.
doPar = logical(sa.useParallel);
if doPar
    try
        doPar = license('test','Distrib_Computing_Toolbox') ~= 0;
    catch
        doPar = false;
    end
end

poolStartedHere = false;
if doPar
    try
        ppool = gcp('nocreate');
        if ~isempty(ppool)
            profNow = string(ppool.Cluster.Profile);
            if ~strcmpi(profNow, 'Processes')
                delete(ppool);
                ppool = [];
            end
        end
        if isempty(ppool)
            if isempty(sa.numWorkers)
                ppool = parpool('Processes');
            else
                ppool = parpool('Processes', sa.numWorkers);
            end
            poolStartedHere = true;
        end
        fprintf('Parallel pool active: profile=%s | workers=%d\n', string(ppool.Cluster.Profile), ppool.NumWorkers);
    catch ME
        fprintf('Could not start parallel pool (%s). Falling back to serial.\n', ME.message);
        doPar = false;
    end
end

% Master seed.
if isempty(sa.seedInit)
    rng('shuffle');
    masterSeed = randi([0, 2^32-1], 1, 1, 'uint32');
else
    masterSeed = uint32(sa.seedInit);
end
fprintf('Weighted-GS RNG masterSeed (Threefry) = %u\n', masterSeed);

% Initial x candidates.
sX0 = RandStream('Threefry','Seed',double(masterSeed));
X0 = cell(nStarts,1);
X0{1} = x_baseline(:);
for k = 2:nStarts
    sX0.Substream = k;
    X0{k} = x_baseline(:) + 0.2 * randn(sX0, numel(x_baseline), 1);
end

results = repmat(struct( ...
    'startIndex', NaN, ...
    'bestMI', -Inf, ...
    'x_best', [], ...
    'runtime', NaN, ...
    'masterSeed', masterSeed, ...
    'substream', NaN), nStarts, 1);

if doPar
    constantStream = parallel.pool.Constant(@() RandStream('Threefry','Seed',double(masterSeed)));
else
    constantStream = [];
end

p_fixed = p_fixed(:);

if doPar
    parfor k = 1:nStarts
        label = sprintf('%s-%02d', taskLabel, k);
        s = constantStream.Value;
        s.Substream = k;
        RandStream.setGlobalStream(s);

        t0 = tic;
        outk = simulated_annealing_gs_weighted(cfg, X0{k}, p_fixed, k, label);
        rt = toc(t0);

        results(k).startIndex = k;
        results(k).bestMI = outk.bestMI;
        results(k).x_best = outk.x_best(:);
        results(k).runtime = rt;
        results(k).substream = k;
    end
else
    s = RandStream('Threefry','Seed',double(masterSeed));
    RandStream.setGlobalStream(s);
    for k = 1:nStarts
        label = sprintf('%s-%02d', taskLabel, k);
        s.Substream = k;
        RandStream.setGlobalStream(s);

        t0 = tic;
        outk = simulated_annealing_gs_weighted(cfg, X0{k}, p_fixed, k, label);
        rt = toc(t0);

        results(k).startIndex = k;
        results(k).bestMI = outk.bestMI;
        results(k).x_best = outk.x_best(:);
        results(k).runtime = rt;
        results(k).substream = k;
    end
end

allMI = [results.bestMI];
[bestMI, idx] = max(allMI);

out = struct();
out.bestStart = idx;
out.bestMI = bestMI;
out.bestX = results(idx).x_best(:);
out.masterSeed = masterSeed;

if doPar && poolStartedHere && isfield(sa,'closePoolWhenDone') && sa.closePoolWhenDone
    delete(gcp('nocreate'));
end
end
