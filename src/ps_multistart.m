function [out, results] = ps_multistart(cfg, x_fixed, p_baseline, taskLabel)
%PS_MULTISTART  Multi-start SA for probabilistic shaping.
%
%   [out, results] = ps_multistart(cfg, x_fixed, p_baseline, taskLabel)
%
%   Similar in spirit to sa_multistart.m, but optimizes probabilities p for
%   a fixed constellation x_fixed.
%
%   Required:
%       cfg.PS.nStarts
%       cfg.PS.useParallel
%       cfg.PS.seedInit
%       cfg.PS.p_min
%       cfg.PS_Evaluator

    if nargin < 4 || isempty(taskLabel)
        taskLabel = '[PS]';
    end
    if isstring(taskLabel), taskLabel = char(taskLabel); end

    ps = cfg.PS;
    x_fixed = x_fixed(:);
    p_baseline = project_probabilities_power(p_baseline, x_fixed, cfg.P_avg, ps.p_min);

    if ~isfield(ps,'nStarts') || isempty(ps.nStarts), ps.nStarts = 1; end
    if ~isfield(ps,'useParallel') || isempty(ps.useParallel), ps.useParallel = false; end
    if ~isfield(ps,'seedInit'), ps.seedInit = []; end
    if ~isfield(ps,'numWorkers'), ps.numWorkers = []; end
    if ~isfield(ps,'closePoolWhenDone'), ps.closePoolWhenDone = false; end

    nStarts = ps.nStarts;
    fprintf('Running PS multistart: nStarts=%d (parallel=%d)\n', nStarts, ps.useParallel);

    doPar = logical(ps.useParallel);
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
                if isempty(ps.numWorkers)
                    ppool = parpool('Processes');
                else
                    ppool = parpool('Processes', ps.numWorkers);
                end
                poolStartedHere = true;
            end

            fprintf('Parallel pool active: profile=%s | workers=%d\n', ...
                string(ppool.Cluster.Profile), ppool.NumWorkers);
        catch ME
            fprintf('Could not start parallel pool (%s). Falling back to serial.\n', ME.message);
            doPar = false;
        end
    end

    if isempty(ps.seedInit)
        rng('shuffle');
        masterSeed = randi([0, 2^32-1], 1, 1, 'uint32');
    else
        masterSeed = uint32(ps.seedInit);
    end
    fprintf('PS RNG masterSeed (Threefry) = %u\n', masterSeed);

    sP0 = RandStream('Threefry','Seed',double(masterSeed));
    P0 = cell(nStarts,1);
    P0{1} = p_baseline(:);

    noiseStd = getfield_default(ps, 'initNoiseStd', 0.05);
    for k = 2:nStarts
        sP0.Substream = k;
        p0k = p_baseline(:) + noiseStd * randn(sP0, numel(p_baseline), 1);
        P0{k} = project_probabilities_power(p0k, x_fixed, cfg.P_avg, ps.p_min);
    end

    results = repmat(struct( ...
        'startIndex', NaN, ...
        'bestMI', -Inf, ...
        'p_best', [], ...
        'runtime', NaN, ...
        'masterSeed', masterSeed, ...
        'substream', NaN), nStarts, 1);

    if doPar
        constantStream = parallel.pool.Constant(@() RandStream('Threefry','Seed',double(masterSeed)));

        parfor k = 1:nStarts
            label = sprintf('%s-%02d', taskLabel, k);

            s = constantStream.Value;
            s.Substream = k;
            RandStream.setGlobalStream(s);

            t0 = tic;
            outk = simulated_annealing_ps(cfg, x_fixed, P0{k}, k, label);
            rt = toc(t0);

            results(k).startIndex = k;
            results(k).bestMI     = outk.bestMI;
            results(k).p_best     = outk.p_best(:);
            results(k).runtime    = rt;
            results(k).substream  = k;
        end
    else
        s = RandStream('Threefry','Seed',double(masterSeed));
        RandStream.setGlobalStream(s);

        for k = 1:nStarts
            label = sprintf('%s-%02d', taskLabel, k);

            s.Substream = k;
            RandStream.setGlobalStream(s);

            t0 = tic;
            outk = simulated_annealing_ps(cfg, x_fixed, P0{k}, k, label);
            rt = toc(t0);

            results(k).startIndex = k;
            results(k).bestMI     = outk.bestMI;
            results(k).p_best     = outk.p_best(:);
            results(k).runtime    = rt;
            results(k).substream  = k;
        end
    end

    allMI = [results.bestMI];
    [bestMI, idx] = max(allMI);

    out = struct();
    out.bestStart  = idx;
    out.bestMI     = bestMI;
    out.bestP      = results(idx).p_best(:);
    out.masterSeed = masterSeed;

    if doPar && poolStartedHere && ps.closePoolWhenDone
        delete(gcp('nocreate'));
    end
end


function val = getfield_default(s, field, defaultVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end
