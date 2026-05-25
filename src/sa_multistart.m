function [out, results] = sa_multistart(cfg, x_baseline, taskLabel)
%SA_MULTISTART  Multi-start SA with optional PARFOR parallelism.
%
% - Forces PROCESS-based parallel pool only.
% - Uses a single master seed (random by default) + substreams per start.
% - Prints the master seed for reproducibility/debug.
%


if nargin < 3 || isempty(taskLabel)
    taskLabel = "[SA]";
end
if isstring(taskLabel), taskLabel = char(taskLabel); end

sa = cfg.SA;

% ---- Defaults (safe) ----
if ~isfield(sa,'nStarts') || isempty(sa.nStarts), sa.nStarts = 1; end
if ~isfield(sa,'useParallel') || isempty(sa.useParallel), sa.useParallel = false; end

% If seedInit is empty -> random master seed each run (client-side).
% If seedInit is numeric -> reproducible run.
if ~isfield(sa,'seedInit'), sa.seedInit = []; end

% Worker count: [] => default
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

        % If a pool exists but is not Processes, restart as Processes
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

        fprintf('Parallel pool active: profile=%s | workers=%d\n', string(p.Cluster.Profile), p.NumWorkers);
    catch ME
        fprintf('Could not start parallel pool (%s). Falling back to serial.\n', ME.message);
        doPar = false;
        p = [];
    end
end

% ---- Choose master seed (client-side only) ----
if isempty(sa.seedInit)
    % random each run, but done ONLY on client (safe)
    rng("shuffle");
    masterSeed = randi([0, 2^32-1], 1, 1, "uint32");
else
    % reproducible mode
    masterSeed = uint32(sa.seedInit);
end

fprintf('RNG masterSeed (Threefry) = %u\n', masterSeed);

% ---- Build initial points X0 using the same masterSeed (reproducible) ----
% Use a separate stream so SA random stream state does not affect X0 generation.
sX0 = RandStream('Threefry','Seed',double(masterSeed));
X0 = cell(nStarts,1);
X0{1} = x_baseline(:);
for k = 2:nStarts
    sX0.Substream = k;   % deterministic per start
    X0{k} = x_baseline(:) + 0.2*randn(sX0, numel(x_baseline), 1);
end

% ---- Preallocate results ----
results = repmat(struct( ...
    'startIndex', NaN, ...
    'bestMI', -Inf, ...
    'x_best', [], ...
    'runtime', NaN, ...
    'masterSeed', masterSeed, ...
    'substream', NaN), nStarts, 1);

% ---- Constant stream for SA randomness (one per worker, lightweight) ----
if doPar
    constantStream = parallel.pool.Constant(@() RandStream('Threefry','Seed',double(masterSeed)));
else
    constantStream = []; 
end

% ---- Run starts ----
if doPar
    parfor k = 1:nStarts
        label = sprintf("%s-%02d", taskLabel, k);

        % Each iteration uses its own substream k (independent of worker scheduling)
        s = constantStream.Value;
        s.Substream = k;
        RandStream.setGlobalStream(s);

        t0 = tic;
        outk = simulated_annealing(cfg, X0{k}, k, label);
        rt = toc(t0);

        results(k).startIndex = k;
        results(k).bestMI     = outk.bestMI;
        results(k).x_best     = outk.x_best(:);
        results(k).runtime    = rt;
        results(k).substream  = k;
    end
else
    % Serial: still use the same scheme for consistency
    s = RandStream('Threefry','Seed',double(masterSeed));
    RandStream.setGlobalStream(s);

    for k = 1:nStarts
        label = sprintf("%s-%02d", taskLabel, k);

        s.Substream = k;
        RandStream.setGlobalStream(s);

        t0 = tic;
        outk = simulated_annealing(cfg, X0{k}, k, label);
        rt = toc(t0);

        results(k).startIndex = k;
        results(k).bestMI     = outk.bestMI;
        results(k).x_best     = outk.x_best(:);
        results(k).runtime    = rt;
        results(k).substream  = k;
    end
end

% ---- Pick best ----
allMI = [results.bestMI];
[bestMI, idx] = max(allMI);

out = struct();
out.bestStart  = idx;
out.bestMI     = bestMI;
out.bestX      = results(idx).x_best(:);
out.masterSeed = masterSeed;

% ---- Optional pool cleanup ----
if doPar && poolStartedHere && isfield(sa,'closePoolWhenDone') && sa.closePoolWhenDone
    delete(gcp("nocreate"));
end

end
