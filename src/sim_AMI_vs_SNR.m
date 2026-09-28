function bundle = sim_AMI_vs_SNR(M, varargin)
%SIM_AMI_VS_SNR  Reproducible AMI-vs-SNR sweep for the IEEE revision.
%
%   bundle = sim_AMI_vs_SNR(M)
%   bundle = sim_AMI_vs_SNR(M, 'Name', Value, ...)
%
% The driver uses the canonical build_fso_config(), validated multi-start
% selection, deterministic per-case seeds, per-case result files, and one
% lightweight summary bundle consumed by plot_AMI_vs_SNR_allM.m.
%
% Name-value options:
%   'P_avg'           average optical intensity (default 1)
%   'SNRVec'          SNR grid in dB (default 5:5:30)
%   'TurbulenceVec'   normalized intensity variance (default [0 .1 .2 .3])
%   'ghN_h'           GH order (default 40)
%   'saMaxIter'       SA iterations per restart (default 10000)
%   'saNStarts'       number of SA restarts (default 8)
%   'minGap'          minimum adjacent spacing (default 0.05)
%   'pinZero'         force first intensity level to zero (default true)
%   'baseSeed'        deterministic experiment seed (default 20260928)
%   'historyEvery'    convergence-history cadence ([] -> temperature block)
%   'UseParallelGrid' parallelize channel cases (default true)
%   'ComputeBER'      auxiliary BER diagnostic (default false)
%   'ResultsRoot'     output root (default repository/results/ieee_revision)
%
% IMPORTANT: SA restarts are serial inside each outer PARFOR case to avoid
% nested parallelism. Each restart winner is validated by sa_multistart().

    if nargin < 1 || isempty(M), M = 4; end

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'M', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 2 && mod(v,1)==0);
    addParameter(p, 'P_avg', 1, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    addParameter(p, 'SNRVec', 5:5:30, @(v) isnumeric(v) && isvector(v) && all(isfinite(v)));
    addParameter(p, 'TurbulenceVec', [0 0.1 0.2 0.3], @(v) isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>=0));
    addParameter(p, 'ghN_h', 40, @(v) isnumeric(v) && isscalar(v) && v>=2 && mod(v,1)==0);
    addParameter(p, 'saMaxIter', 10000, @(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p, 'saNStarts', 8, @(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p, 'minGap', 0.05, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p, 'pinZero', true, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'baseSeed', 20260928, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p, 'historyEvery', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0));
    addParameter(p, 'UseParallelGrid', true, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'ComputeBER', false, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'ResultsRoot', fso_result_utils.default_results_root(), @(v) ischar(v) || isstring(v));
    parse(p, M, varargin{:});
    o = p.Results;

    M       = double(M);
    P_avg   = double(o.P_avg);
    snrVec  = double(o.SNRVec(:).');
    sigVec  = double(o.TurbulenceVec(:).');
    nSNR    = numel(snrVec);
    nSig    = numel(sigVec);
    nPoints = nSNR * nSig;

    xMaxBound = fso_result_utils.feasible_xmax_bound(M, P_avg, o.minGap, o.pinZero);

    resultsRoot = char(o.ResultsRoot);
    outDir = fullfile(resultsRoot, 'ami_vs_snr', sprintf('M%d', M));
    caseDir = fullfile(outDir, 'cases');
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end

    runMeta = fso_result_utils.run_metadata(mfilename);

    fprintf('\n============================================================\n');
    fprintf('IEEE revision AMI-vs-SNR sweep\n');
    fprintf('M=%d | Pavg=%.4g | minGap=%.4g | xMaxBound=%.4g\n', M, P_avg, o.minGap, xMaxBound);
    fprintf('SNR points=%d | turbulence points=%d | total cases=%d\n', nSNR, nSig, nPoints);
    fprintf('SA: %d restarts x %d iterations | GH=%d\n', o.saNStarts, o.saMaxIter, o.ghN_h);
    fprintf('Validated winner selection enabled for every case.\n');
    fprintf('Results: %s\n', outDir);
    fprintf('============================================================\n\n');

    [SNRGrid, SigGrid] = meshgrid(snrVec, sigVec);
    taskSNR = SNRGrid(:);
    taskSig = SigGrid(:);

    caseResults = cell(nPoints,1);
    caseFiles   = cell(nPoints,1);

    doPar = logical(o.UseParallelGrid);
    if doPar
        try
            doPar = license('test','Distrib_Computing_Toolbox') ~= 0;
        catch
            doPar = false;
        end
    end

    if doPar
        pool = gcp('nocreate');
        if isempty(pool), pool = parpool('Processes'); end
        fprintf('Outer grid parallelism: %d process workers\n', pool.NumWorkers);

        Q = parallel.pool.DataQueue;
        doneCount = 0;
        afterEach(Q, @(~) progress_tick());

        parfor k = 1:nPoints
            [caseResults{k}, caseFiles{k}] = run_case( ...
                M, P_avg, taskSNR(k), taskSig(k), o, xMaxBound, runMeta, caseDir);
            send(Q,1);
        end
    else
        fprintf('Outer grid parallelism: disabled\n');
        for k = 1:nPoints
            [caseResults{k}, caseFiles{k}] = run_case( ...
                M, P_avg, taskSNR(k), taskSig(k), o, xMaxBound, runMeta, caseDir);
            fprintf('  completed %d/%d\n', k, nPoints);
        end
    end

    caseGrid = reshape(caseResults, nSig, nSNR);
    fileGrid = reshape(caseFiles, nSig, nSNR);

    summary = init_summary(nSig, nSNR, logical(o.ComputeBER));
    for iSig = 1:nSig
        for iSNR = 1:nSNR
            c = caseGrid{iSig,iSNR};
            summary.amiPAMFast(iSig,iSNR)      = c.baseline.amiFast;
            summary.amiPAMValidated(iSig,iSNR) = c.baseline.amiValidated;
            summary.amiGSFast(iSig,iSNR)       = c.winner.amiFast;
            summary.amiGSValidated(iSig,iSNR)  = c.winner.amiValidated;
            summary.gainBits(iSig,iSNR)        = c.metrics.gainBits;
            summary.selectionChanged(iSig,iSNR)= c.selection.changed;
            summary.fastWinnerIndex(iSig,iSNR) = c.selection.fastWinnerIndex;
            summary.validatedWinnerIndex(iSig,iSNR) = c.selection.validatedWinnerIndex;
            summary.meanAbsFastValidationGap(iSig,iSNR) = c.selection.meanAbsFastValidationGap;
            summary.maxAbsFastValidationGap(iSig,iSNR)  = c.selection.maxAbsFastValidationGap;
            summary.saRuntime(iSig,iSNR)        = c.runtime.sa;
            summary.validationRuntime(iSig,iSNR)= c.runtime.validation;
            summary.xGS{iSig,iSNR}              = c.winner.x;
            if o.ComputeBER
                summary.berPAM(iSig,iSNR) = c.baseline.ber;
                summary.berGS(iSig,iSNR)  = c.winner.ber;
            end
        end
    end

    bundle = struct();
    bundle.meta       = runMeta;
    bundle.experiment = 'AMI_vs_SNR';
    bundle.axes       = struct('M',M, 'P_avg',P_avg, 'SNR_dB',snrVec, 'sigma_X_sq',sigVec);
    bundle.settings   = struct( ...
        'ghN_h',o.ghN_h, 'saMaxIter',o.saMaxIter, 'saNStarts',o.saNStarts, ...
        'minGap',o.minGap, 'pinZero',o.pinZero, 'baseSeed',o.baseSeed, ...
        'historyEvery',o.historyEvery, 'computeBER',o.ComputeBER, ...
        'xMaxBound',xMaxBound);
    bundle.caseFiles = fileGrid;
    bundle.summary   = summary;

    bundlePath = fullfile(outDir, sprintf('AMI_vs_SNR_M%d.mat', M));
    save(bundlePath, 'bundle', '-v7.3');

    fprintf('\nSaved summary bundle: %s\n', bundlePath);
    print_summary(bundle);

    function progress_tick()
        doneCount = doneCount + 1;
        fprintf('  completed %d/%d\n', doneCount, nPoints);
    end
end


function [caseResult, casePath] = run_case(M, P_avg, SNR_dB, sigma_X_sq, o, xMaxBound, runMeta, caseDir)
    seed = fso_result_utils.case_seed(o.baseSeed, M, SNR_dB, sigma_X_sq);

    cfg = build_fso_config(M, P_avg, SNR_dB, sigma_X_sq, ...
        'ghN_h', o.ghN_h, ...
        'xMaxBound', xMaxBound, ...
        'saMaxIter', o.saMaxIter, ...
        'saNStarts', o.saNStarts, ...
        'saUseParallel', false, ...
        'minGap', o.minGap, ...
        'pinZero', o.pinZero, ...
        'seedInit', seed, ...
        'logEvery', 0, ...
        'historyEvery', o.historyEvery);

    xPAM = define_constellation(M, P_avg, true, "mean");
    AMI_functions.assert_constellation_feasible(xPAM, cfg, 1e-10);

    amiPAMFast = cfg.AMI_Evaluator(xPAM);
    tBaseVal = tic;
    amiPAMVal = cfg.AMI_Validator(xPAM);
    baseValidationRuntime = toc(tBaseVal);

    label = sprintf('[M%d|SNR%.1f|sig%.3f]', M, SNR_dB, sigma_X_sq);
    [saOut, starts] = sa_multistart(cfg, xPAM, label);

    xGS = saOut.bestXValidated(:);
    amiGSVal  = saOut.bestMIValidated;
    amiGSFast = saOut.validatedWinnerFastMI;

    berPAM = NaN;
    berGS  = NaN;
    if o.ComputeBER
        berPAM = BER_functions.calculate_BER_noCSI_ML(xPAM(:), cfg);
        berGS  = BER_functions.calculate_BER_noCSI_ML(xGS(:), cfg);
    end

    saOutSaved = saOut;
    if isfield(saOutSaved,'bestRun'), saOutSaved = rmfield(saOutSaved,'bestRun'); end

    caseResult = struct();
    caseResult.meta = runMeta;
    caseResult.case = struct('M',M, 'P_avg',P_avg, 'SNR_dB',SNR_dB, ...
        'sigma_X_sq',sigma_X_sq, 'seed',seed);
    caseResult.config = fso_result_utils.config_snapshot(cfg);
    caseResult.baseline = struct('x',xPAM(:), 'amiFast',amiPAMFast, ...
        'amiValidated',amiPAMVal, 'ber',berPAM);
    caseResult.starts = starts;
    caseResult.optimizer = saOutSaved;
    caseResult.winner = struct( ...
        'startIndex',saOut.bestStartValidated, ...
        'x',xGS, ...
        'amiFast',amiGSFast, ...
        'amiValidated',amiGSVal, ...
        'ber',berGS);
    caseResult.selection = struct( ...
        'fastWinnerIndex',saOut.bestStartFast, ...
        'validatedWinnerIndex',saOut.bestStartValidated, ...
        'changed',saOut.selectionChanged, ...
        'validatedAdvantageOverFastWinner',saOut.validatedAdvantageOverFastWinner, ...
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap, ...
        'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    caseResult.metrics = struct( ...
        'gainBits',amiGSVal-amiPAMVal, ...
        'fastGainBits',amiGSFast-amiPAMFast);
    caseResult.runtime = struct( ...
        'sa',saOut.totalSARuntime, ...
        'validation',saOut.validationTotalRuntime + baseValidationRuntime, ...
        'total',saOut.totalRuntime + baseValidationRuntime);

    casePath = fullfile(caseDir, fso_result_utils.case_filename(M,SNR_dB,sigma_X_sq));
    save_case(casePath, caseResult);
end


function save_case(path, caseResult)
% Separate helper keeps PARFOR file output transparent and case-specific.
    save(path, 'caseResult', '-v7.3');
end


function s = init_summary(nSig, nSNR, computeBER)
    s = struct();
    s.amiPAMFast      = nan(nSig,nSNR);
    s.amiPAMValidated = nan(nSig,nSNR);
    s.amiGSFast       = nan(nSig,nSNR);
    s.amiGSValidated  = nan(nSig,nSNR);
    s.gainBits        = nan(nSig,nSNR);
    s.selectionChanged = false(nSig,nSNR);
    s.fastWinnerIndex = nan(nSig,nSNR);
    s.validatedWinnerIndex = nan(nSig,nSNR);
    s.meanAbsFastValidationGap = nan(nSig,nSNR);
    s.maxAbsFastValidationGap  = nan(nSig,nSNR);
    s.saRuntime         = nan(nSig,nSNR);
    s.validationRuntime = nan(nSig,nSNR);
    s.xGS = cell(nSig,nSNR);
    if computeBER
        s.berPAM = nan(nSig,nSNR);
        s.berGS  = nan(nSig,nSNR);
    end
end


function print_summary(bundle)
    snrVec = bundle.axes.SNR_dB;
    sigVec = bundle.axes.sigma_X_sq;
    S = bundle.summary;

    fprintf('\nValidated AMI shaping gain [bits/symbol]\n');
    fprintf('sigma_X^2 \\ SNR');
    for j = 1:numel(snrVec), fprintf(' | %6.1f', snrVec(j)); end
    fprintf('\n');
    for i = 1:numel(sigVec)
        fprintf('%12.3f', sigVec(i));
        for j = 1:numel(snrVec), fprintf(' | %+6.4f', S.gainBits(i,j)); end
        fprintf('\n');
    end

    nChanged = nnz(S.selectionChanged);
    fprintf('\nFast-vs-validated winner changed in %d/%d cases.\n', ...
        nChanged, numel(S.selectionChanged));
    fprintf('Maximum |fast - validated| among selected candidates: %.3e bits/symbol\n', ...
        max(S.maxAbsFastValidationGap(:), [], 'omitnan'));
end
