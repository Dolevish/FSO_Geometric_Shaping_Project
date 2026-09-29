function bundle = sim_fast_evaluator_optimizer_benchmark(source,varargin)
%SIM_FAST_EVALUATOR_OPTIMIZER_BENCHMARK  Commit-I.4 GH400 vs log-h benchmark.
%
% This experiment compares the historical Commit-G Gauss-Hermite fast AMI
% objective against the Commit-I.3 production uniform-log-h evaluator on
% exactly the same representative FSO case and, by default, the same SA
% budget/seeds used by the supplied Commit-G result bundle.
%
% The experiment has two layers:
%   A) evaluator-only: rescore EVERY restart winner from the reference GH run
%      with GH, log-h, and the independent validator, including repeated
%      runtime measurements on the exact same constellations/y-grid;
%   B) optimizer-level: rerun SA with log-h only, reusing the historical
%      replicate master seeds and the selected Commit-G SA profile.
%
% SOURCE must be the path to saHyperparameterSensitivity.mat from the
% completed Commit-G GH reference run. The reference case files are loaded
% from the source bundle or, if paths moved between computers, from the
% sibling cases/ directory.
%
% Default scientific target:
%   M=32, SNR=20 dB, sigma_X^2=0.1, d_min=0.01
%   historical reference: GH400
%   log-h: dt=0.01, span=mu_t +/- 8 sigma_t, yBlockSize=512
%   SA budget/profile: inherited from SOURCE (normally 4 x 4000, 2 reps,
%   T0=0.4, Tf=1e-3, step0=0.15, block=50).
%
% Name-value options:
%   'ReferenceGHOrder'  expected GH order in source (default 400)
%   'loghDt'            log-h spacing (default 0.01)
%   'loghSpanSigma'     half-span in sigma_t (default 8)
%   'loghYBlockSize'    y block size (default 512)
%   'TimingWarmups'     warm-up evaluations per method/candidate (default 1)
%   'TimingRepeats'     timed evaluations per method/candidate (default 3)
%   'RunSA'             run end-to-end log-h SA layer (default true)
%   'ReplicateVec'      source replicate indices ([] -> all)
%   'saMaxIter'         [] -> inherit source; override only for smoke tests
%   'saNStarts'         [] -> inherit source; override only for smoke tests
%   'UseParallelCases'  parallelize log-h replicate cases (default true)
%   'MakePlots'         create figures (default true)
%   'SaveFigures'       save PNG/FIG pairs (default true)
%   'RunLabel'          optional run label
%   'ResultsRoot'       output root
%
% A budget override is intentionally marked as NON-COMPARABLE to the stored
% GH optimizer result. It is useful only for smoke testing the code path.

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p,'source',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ReferenceGHOrder',400,@positive_integer);
    addParameter(p,'loghDt',0.01,@positive_scalar);
    addParameter(p,'loghSpanSigma',8,@positive_scalar);
    addParameter(p,'loghYBlockSize',512,@positive_integer);
    addParameter(p,'TimingWarmups',1,@nonnegative_integer);
    addParameter(p,'TimingRepeats',3,@positive_integer);
    addParameter(p,'RunSA',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'ReplicateVec',[],@(v) isempty(v) || positive_integer_vector(v));
    addParameter(p,'saMaxIter',[],@(v) isempty(v) || positive_integer(v));
    addParameter(p,'saNStarts',[],@(v) isempty(v) || positive_integer(v));
    addParameter(p,'UseParallelCases',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'MakePlots',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'SaveFigures',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v) || isstring(v));
    parse(p,source,varargin{:});
    o = p.Results;

    [src,sourcePath,sourceDir] = load_reference_bundle(source);
    validate_reference_bundle(src);

    baselineIdx = double(src.summary.baselineDefinitionIndex);
    sourceNRep = double(src.settings.nReplicates);
    sourceNStarts = double(src.settings.saNStarts);
    sourceMaxIter = double(src.settings.saMaxIter);

    if isempty(o.ReplicateVec)
        repVec = 1:sourceNRep;
    else
        repVec = unique(double(o.ReplicateVec(:).'),'stable');
        if any(repVec > sourceNRep)
            error('sim_fast_evaluator_optimizer_benchmark:BadReplicate', ...
                'ReplicateVec exceeds source nReplicates=%d.',sourceNRep);
        end
    end
    nRep = numel(repVec);

    M = double(src.case.M);
    Pavg = double(src.case.P_avg);
    SNR_dB = double(src.case.SNR_dB);
    sigma_X_sq = double(src.case.sigma_X_sq);
    minGap = double(src.case.minGap);
    pinZero = true; % Commit-G experiment hard-coded pinZero=true.

    if ~isfield(src.case,'ghN_h') || double(src.case.ghN_h) ~= double(o.ReferenceGHOrder)
        error('sim_fast_evaluator_optimizer_benchmark:ReferenceGHMismatch', ...
            'Source GH order does not match requested ReferenceGHOrder=%d.',o.ReferenceGHOrder);
    end

    if isfield(src.case,'fastFadingMethod') && ~strcmpi(src.case.fastFadingMethod,'gh')
        error('sim_fast_evaluator_optimizer_benchmark:ReferenceNotGH', ...
            'Source explicitly reports fastFadingMethod=%s, expected gh.',src.case.fastFadingMethod);
    end

    saMaxIter = sourceMaxIter;
    saNStarts = sourceNStarts;
    if ~isempty(o.saMaxIter), saMaxIter = double(o.saMaxIter); end
    if ~isempty(o.saNStarts), saNStarts = double(o.saNStarts); end
    budgetComparable = (saMaxIter == sourceMaxIter) && (saNStarts == sourceNStarts);

    historyEvery = double(src.settings.historyEvery);
    saProfile = src.baselineSA;
    xPAM = real(src.baselinePAM.x(:));
    xMaxBound = fso_result_utils.feasible_xmax_bound(M,Pavg,minGap,pinZero);

    refCases = cell(nRep,1);
    refCasePaths = cell(nRep,1);
    for r = 1:nRep
        [refCases{r},refCasePaths{r}] = load_reference_case( ...
            src,sourceDir,baselineIdx,repVec(r));
        validate_reference_case(refCases{r},repVec(r),M,SNR_dB,sigma_X_sq,minGap);
        if numel(refCases{r}.starts) ~= sourceNStarts
            error('sim_fast_evaluator_optimizer_benchmark:StartCountMismatch', ...
                'Reference replicate %d contains %d starts, expected %d.', ...
                repVec(r),numel(refCases{r}.starts),sourceNStarts);
        end
    end

    commonCfgArgs = { ...
        'xMaxBound',xMaxBound,'saMaxIter',1,'saNStarts',1, ...
        'saUseParallel',false,'minGap',minGap,'pinZero',pinZero, ...
        'seedInit',1,'logEvery',0,'historyEvery',1};

    cfgGH = build_fso_config(M,Pavg,SNR_dB,sigma_X_sq,commonCfgArgs{:}, ...
        'FastFadingMethod','gh','ghN_h',o.ReferenceGHOrder);
    cfgLoghEval = build_fso_config(M,Pavg,SNR_dB,sigma_X_sq,commonCfgArgs{:}, ...
        'FastFadingMethod','logh','loghDt',o.loghDt, ...
        'loghSpanSigma',o.loghSpanSigma,'loghYBlockSize',o.loghYBlockSize);

    if numel(cfgGH.y_grid) ~= numel(cfgLoghEval.y_grid) || ...
            max(abs(cfgGH.y_grid(:)-cfgLoghEval.y_grid(:))) > 1e-14
        error('sim_fast_evaluator_optimizer_benchmark:YGridMismatch', ...
            'GH and log-h configs do not use the same outer y-grid.');
    end

    AMI_functions.assert_constellation_feasible(xPAM,cfgGH,1e-10);
    baselineValidatorNow = cfgGH.AMI_Validator(xPAM);
    baselineGHNow = cfgGH.AMI_Evaluator(xPAM);
    baselineLoghNow = cfgLoghEval.AMI_Evaluator(xPAM);

    % ------------------------------------------------------------------
    % Immutable result directory.
    % ------------------------------------------------------------------
    conditionTag = sprintf('M%d_SNR%s_sig%s_gap%s', ...
        M,num_token(SNR_dB,2),num_token(sigma_X_sq,4),num_token(minGap,4));
    parentDir = fullfile(char(o.ResultsRoot),'fast_evaluator_optimizer_benchmark',conditionTag);
    labelToken = sanitize_run_label(o.RunLabel);
    budgetToken = sprintf('LOGHdt%s_GH%d_I%d_S%d_R%d', ...
        num_token(o.loghDt,4),o.ReferenceGHOrder,saMaxIter,saNStarts,nRep);
    utcToken = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requestedRunId = sprintf('%s_%s_%s',labelToken,budgetToken,utcToken);
    [outDir,runId] = make_unique_run_directory(parentDir,requestedRunId);
    caseDir = fullfile(outDir,'cases'); mkdir(caseDir);
    figDir = fullfile(outDir,'figures');
    if o.MakePlots && o.SaveFigures, mkdir(figDir); end

    runMeta = fso_result_utils.run_metadata(mfilename);
    runMeta.runId = runId;
    runMeta.runLabel = char(string(o.RunLabel));
    runMeta.outputDirectory = outDir;
    runMeta.sourcePath = sourcePath;

    fprintf('\n============================================================\n');
    fprintf('Commit I.4 - GH%d vs production log-h optimizer benchmark\n',o.ReferenceGHOrder);
    fprintf('Run ID: %s\n',runId);
    fprintf('Source: %s\n',sourcePath);
    fprintf('Case: M=%d | SNR=%.2f dB | sigma_X^2=%.4f | d_min=%.4f\n', ...
        M,SNR_dB,sigma_X_sq,minGap);
    fprintf('Reference evaluator: GH%d | log-h: dt=%.5g, span=+/-%.3g sigma_t, yBlock=%d\n', ...
        o.ReferenceGHOrder,o.loghDt,o.loghSpanSigma,o.loghYBlockSize);
    fprintf('Reference source SA: %d starts x %d iterations | selected reps=%s\n', ...
        sourceNStarts,sourceMaxIter,mat2str(repVec));
    fprintf('Requested log-h SA: %d starts x %d iterations | budget comparable=%d\n', ...
        saNStarts,saMaxIter,budgetComparable);
    fprintf('Evaluator timing: warmups=%d | repeats=%d\n',o.TimingWarmups,o.TimingRepeats);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    % ------------------------------------------------------------------
    % Layer A: evaluator-only benchmark on all historical restart winners.
    % ------------------------------------------------------------------
    evalBench = evaluator_benchmark(refCases,repVec,cfgGH,cfgLoghEval, ...
        double(o.TimingWarmups),double(o.TimingRepeats),sourceNStarts);

    fprintf('Evaluator-only benchmark complete: %d constellations.\n',numel(evalBench.validatorAMI));
    fprintf('  mean |GH%d-validator|  = %.3e\n',o.ReferenceGHOrder,mean(abs(evalBench.ghError(:)),'omitnan'));
    fprintf('  mean |logh-validator|   = %.3e\n',mean(abs(evalBench.loghError(:)),'omitnan'));
    fprintf('  median GH%d time/eval   = %.4fs\n',o.ReferenceGHOrder,median(evalBench.ghRuntimeMedian(:),'omitnan'));
    fprintf('  median log-h time/eval  = %.4fs\n',median(evalBench.loghRuntimeMedian(:),'omitnan'));
    fprintf('  median GH/log-h speedup = %.3fx\n',median(evalBench.speedupGHOverLogh(:),'omitnan'));
    fprintf('  max |recomputed GH - stored GH fast| = %.3e\n\n', ...
        max(abs(evalBench.ghReproductionError(:)),[],'omitnan'));

    % ------------------------------------------------------------------
    % Layer B: rerun SA with production log-h only, reusing source seeds.
    % ------------------------------------------------------------------
    loghCases = cell(nRep,1);
    loghCaseFiles = cell(nRep,1);
    doPar = false;

    if o.RunSA
        doPar = logical(o.UseParallelCases);
        if doPar
            try
                doPar = license('test','Distrib_Computing_Toolbox') ~= 0;
            catch
                doPar = false;
            end
        end

        if ~budgetComparable
            fprintf(['WARNING: SA budget override is active. Optimizer-level GH/log-h ' ...
                'performance is a smoke diagnostic, not a controlled A/B comparison.\n']);
        end

        if doPar
            pool = gcp('nocreate');
            if isempty(pool), pool = parpool('Processes'); end
            fprintf('Outer replicate parallelism: %d process workers (tasks=%d)\n', ...
                pool.NumWorkers,nRep);
            Q = parallel.pool.DataQueue;
            doneCount = 0;
            afterEach(Q,@(~) progress_tick());
            parfor r = 1:nRep
                [loghCases{r},loghCaseFiles{r}] = run_logh_case( ...
                    refCases{r},repVec(r),M,Pavg,SNR_dB,sigma_X_sq,minGap,pinZero, ...
                    xMaxBound,xPAM,baselineValidatorNow,saProfile,saMaxIter,saNStarts, ...
                    historyEvery,o,runMeta,caseDir);
                send(Q,1);
            end
        else
            fprintf('Outer replicate parallelism: disabled\n');
            for r = 1:nRep
                [loghCases{r},loghCaseFiles{r}] = run_logh_case( ...
                    refCases{r},repVec(r),M,Pavg,SNR_dB,sigma_X_sq,minGap,pinZero, ...
                    xMaxBound,xPAM,baselineValidatorNow,saProfile,saMaxIter,saNStarts, ...
                    historyEvery,o,runMeta,caseDir);
                fprintf('  completed log-h replicate %d/%d\n',r,nRep);
            end
        end
    else
        fprintf('RunSA=false: optimizer-level layer skipped.\n');
    end

    summary = build_summary(refCases,loghCases,repVec,evalBench, ...
        baselineValidatorNow,budgetComparable,logical(o.RunSA));

    bundle = struct();
    bundle.meta = runMeta;
    bundle.experiment = 'fast_evaluator_optimizer_benchmark';
    bundle.run = struct('id',runId,'label',char(string(o.RunLabel)), ...
        'outputDirectory',outDir,'budgetComparable',budgetComparable);
    bundle.source = struct('path',sourcePath,'runId',source_run_id(src), ...
        'baselineDefinitionIndex',baselineIdx,'caseFiles',{refCasePaths});
    bundle.case = struct('M',M,'P_avg',Pavg,'SNR_dB',SNR_dB, ...
        'sigma_X_sq',sigma_X_sq,'minGap',minGap,'pinZero',pinZero, ...
        'xMaxBound',xMaxBound,'replicateVec',repVec);
    bundle.methods = struct( ...
        'reference',struct('name','gh','order',double(o.ReferenceGHOrder)), ...
        'production',struct('name','logh','dt',double(o.loghDt), ...
            'spanSigma',double(o.loghSpanSigma),'yBlockSize',double(o.loghYBlockSize)), ...
        'sameOuterYGrid',true,'independentValidator',true);
    bundle.sa = struct('sourceMaxIter',sourceMaxIter,'sourceNStarts',sourceNStarts, ...
        'maxIter',saMaxIter,'nStarts',saNStarts,'profile',saProfile, ...
        'historyEvery',historyEvery,'budgetComparable',budgetComparable, ...
        'pairedMasterSeeds',true,'outerParallelism',doPar);
    bundle.baseline = struct('x',xPAM,'storedValidated',double(src.baselinePAM.amiValidated), ...
        'recomputedValidated',baselineValidatorNow,'ghFast',baselineGHNow, ...
        'loghFast',baselineLoghNow);
    bundle.evaluatorBenchmark = evalBench;
    bundle.loghCaseFiles = loghCaseFiles;
    bundle.summary = summary;
    bundle.figureFiles = {};

    if o.MakePlots
        bundle.figureFiles = make_plots(bundle,refCases,loghCases,figDir,logical(o.SaveFigures));
    end

    bundlePath = fullfile(outDir,'evaluatorOptimizerBenchmark.mat');
    save(bundlePath,'bundle','-v7.3');
    fprintf('\nSaved Commit-I.4 benchmark bundle: %s\n',bundlePath);
    print_summary(bundle);

    function progress_tick()
        doneCount = doneCount + 1;
        fprintf('  completed log-h replicate %d/%d\n',doneCount,nRep);
    end
end


function E = evaluator_benchmark(refCases,repVec,cfgGH,cfgLH,nWarm,nTiming,nStarts)
    nRep = numel(refCases);
    shape = [nRep nStarts];
    fields = {'storedFastAMI','storedValidatedAMI','validatorAMI','ghAMI','loghAMI', ...
        'ghError','loghError','ghReproductionError','validatorReproductionError', ...
        'ghRuntimeMedian','loghRuntimeMedian','speedupGHOverLogh'};
    E = struct();
    for k=1:numel(fields), E.(fields{k}) = nan(shape); end
    E.replicateIndex = repVec(:);
    E.restartIndex = 1:nStarts;

    for r=1:nRep
        starts = refCases{r}.starts;
        for s=1:nStarts
            x = real(starts(s).x_best(:));
            AMI_functions.assert_constellation_feasible(x,cfgGH,1e-10);

            val = cfgGH.AMI_Validator(x);
            gh = cfgGH.AMI_Evaluator(x);
            lh = cfgLH.AMI_Evaluator(x);
            tGH = timed_median(cfgGH.AMI_Evaluator,x,nWarm,nTiming);
            tLH = timed_median(cfgLH.AMI_Evaluator,x,nWarm,nTiming);

            E.storedFastAMI(r,s) = starts(s).bestMIFast;
            E.storedValidatedAMI(r,s) = starts(s).bestMIValidated;
            E.validatorAMI(r,s) = val;
            E.ghAMI(r,s) = gh;
            E.loghAMI(r,s) = lh;
            E.ghError(r,s) = gh-val;
            E.loghError(r,s) = lh-val;
            E.ghReproductionError(r,s) = gh-starts(s).bestMIFast;
            E.validatorReproductionError(r,s) = val-starts(s).bestMIValidated;
            E.ghRuntimeMedian(r,s) = tGH;
            E.loghRuntimeMedian(r,s) = tLH;
            E.speedupGHOverLogh(r,s) = tGH/max(tLH,eps);
        end
    end
end


function tMed = timed_median(fun,x,nWarm,nTiming)
    for k=1:nWarm, fun(x); end
    t = nan(1,nTiming);
    for k=1:nTiming
        tt=tic; fun(x); t(k)=toc(tt);
    end
    tMed=median(t,'omitnan');
end


function [c,path] = run_logh_case(refCase,repIdx,M,Pavg,SNR_dB,sigma_X_sq,minGap,pinZero, ...
        xMaxBound,xPAM,baselineValidated,saProfile,saMaxIter,saNStarts,historyEvery,o,runMeta,caseDir)

    seed = double(refCase.case.seed);
    cfg = build_fso_config(M,Pavg,SNR_dB,sigma_X_sq, ...
        'FastFadingMethod','logh','loghDt',o.loghDt, ...
        'loghSpanSigma',o.loghSpanSigma,'loghYBlockSize',o.loghYBlockSize, ...
        'xMaxBound',xMaxBound,'saMaxIter',saMaxIter,'saNStarts',saNStarts, ...
        'saUseParallel',false,'minGap',minGap,'pinZero',pinZero, ...
        'seedInit',seed,'logEvery',0,'historyEvery',historyEvery);

    cfg.SA.T0 = double(saProfile.T0);
    cfg.SA.Tf = double(saProfile.Tf);
    cfg.SA.baseStd0 = double(saProfile.baseStd0);
    cfg.SA.itersPerTemp = double(saProfile.itersPerTemp);
    cfg.SA.nBlocks = ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate = exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);

    label = sprintf('[I4|logh|rep%d]',repIdx);
    tCase=tic;
    [saOut,starts] = sa_multistart(cfg,xPAM,label);
    wall=toc(tCase);

    win = starts(saOut.bestStartValidated);
    k95 = convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.95);
    k99 = convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.99);
    vals = [starts.bestMIValidated];
    nEval = double(saNStarts)*double(saMaxIter);

    c = struct();
    c.meta = runMeta;
    c.case = struct('replicate',repIdx,'seed',seed,'M',M,'SNR_dB',SNR_dB, ...
        'sigma_X_sq',sigma_X_sq,'minGap',minGap);
    c.method = struct('fastFadingMethod','logh','loghDt',double(o.loghDt), ...
        'loghSpanSigma',double(o.loghSpanSigma),'loghYBlockSize',double(o.loghYBlockSize));
    c.sa = struct('T0',cfg.SA.T0,'Tf',cfg.SA.Tf,'baseStd0',cfg.SA.baseStd0, ...
        'itersPerTemp',cfg.SA.itersPerTemp,'coolingRate',cfg.SA.coolingRate, ...
        'maxIter',cfg.SA.maxIter,'nStarts',cfg.SA.nStarts);
    c.config = fso_result_utils.config_snapshot(cfg);
    c.starts = starts;
    c.optimizer = saOut;
    if isfield(c.optimizer,'bestRun'), c.optimizer=rmfield(c.optimizer,'bestRun'); end
    c.winner = struct('startIndex',saOut.bestStartValidated,'x',saOut.bestXValidated(:), ...
        'amiFast',saOut.validatedWinnerFastMI,'amiValidated',saOut.bestMIValidated, ...
        'gainBits',saOut.bestMIValidated-baselineValidated, ...
        'acceptanceRate',win.acceptanceRate,'convergenceIter95',k95, ...
        'convergenceIter99',k99,'history',win.history);
    c.repeatability = struct('restartValidatedMean',mean(vals,'omitnan'), ...
        'restartValidatedStd',std(vals,0,'omitnan'),'restartValidatedRange',max(vals)-min(vals));
    c.validation = struct('selectionChanged',saOut.selectionChanged, ...
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap, ...
        'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    c.runtime = struct('wallSeconds',wall,'saSecondsSummed',saOut.totalSARuntime, ...
        'validationSeconds',saOut.validationTotalRuntime,'fastEvaluations',nEval, ...
        'fastEvaluationsPerSecond',nEval/max(saOut.totalSARuntime,eps));

    fname=sprintf('logh_rep%02d.mat',repIdx);
    path=fullfile(caseDir,fname);
    save(path,'c','-v7.3');
end


function S = build_summary(refCases,loghCases,repVec,E,baselineVal,budgetComparable,runSA)
    nRep=numel(refCases);
    metricNames={'amiValidated','gainBits','conv95','conv99','wallSeconds', ...
        'evalPerSec','acceptanceRate','restartStd','maxValidationGap','selectionChanged'};
    ref=struct(); lh=struct();
    for k=1:numel(metricNames)
        ref.(metricNames{k})=nan(nRep,1);
        lh.(metricNames{k})=nan(nRep,1);
    end
    refX=cell(nRep,1); lhX=cell(nRep,1);

    for r=1:nRep
        c=refCases{r};
        ref.amiValidated(r)=c.winner.amiValidated;
        ref.gainBits(r)=c.winner.amiValidated-baselineVal;
        ref.conv95(r)=c.winner.convergenceIter95;
        ref.conv99(r)=c.winner.convergenceIter99;
        ref.wallSeconds(r)=c.runtime.wallSeconds;
        ref.evalPerSec(r)=c.runtime.fastEvaluationsPerSecond;
        ref.acceptanceRate(r)=c.winner.acceptanceRate;
        ref.restartStd(r)=c.repeatability.restartValidatedStd;
        ref.maxValidationGap(r)=c.validation.maxAbsFastValidationGap;
        ref.selectionChanged(r)=double(c.validation.selectionChanged);
        refX{r}=c.winner.x(:);

        if runSA
            q=loghCases{r};
            lh.amiValidated(r)=q.winner.amiValidated;
            lh.gainBits(r)=q.winner.gainBits;
            lh.conv95(r)=q.winner.convergenceIter95;
            lh.conv99(r)=q.winner.convergenceIter99;
            lh.wallSeconds(r)=q.runtime.wallSeconds;
            lh.evalPerSec(r)=q.runtime.fastEvaluationsPerSecond;
            lh.acceptanceRate(r)=q.winner.acceptanceRate;
            lh.restartStd(r)=q.repeatability.restartValidatedStd;
            lh.maxValidationGap(r)=q.validation.maxAbsFastValidationGap;
            lh.selectionChanged(r)=double(q.validation.selectionChanged);
            lhX{r}=q.winner.x(:);
        end
    end

    S=struct();
    S.replicateVec=repVec;
    S.reference=stats_struct(ref);
    S.logh=stats_struct(lh);
    S.reference.raw=ref;
    S.logh.raw=lh;
    S.reference.xWinner=refX;
    S.logh.xWinner=lhX;
    S.budgetComparable=budgetComparable;

    S.evaluator=struct();
    S.evaluator.meanAbsErrorGH=mean(abs(E.ghError(:)),'omitnan');
    S.evaluator.meanAbsErrorLogh=mean(abs(E.loghError(:)),'omitnan');
    S.evaluator.maxAbsErrorGH=max(abs(E.ghError(:)),[],'omitnan');
    S.evaluator.maxAbsErrorLogh=max(abs(E.loghError(:)),[],'omitnan');
    S.evaluator.meanSignedErrorGH=mean(E.ghError(:),'omitnan');
    S.evaluator.meanSignedErrorLogh=mean(E.loghError(:),'omitnan');
    S.evaluator.medianRuntimeGH=median(E.ghRuntimeMedian(:),'omitnan');
    S.evaluator.medianRuntimeLogh=median(E.loghRuntimeMedian(:),'omitnan');
    S.evaluator.medianSpeedupGHOverLogh=median(E.speedupGHOverLogh(:),'omitnan');
    S.evaluator.maxGHReproductionError=max(abs(E.ghReproductionError(:)),[],'omitnan');
    S.evaluator.maxValidatorReproductionError=max(abs(E.validatorReproductionError(:)),[],'omitnan');

    if runSA
        S.deltaValidatedAMILoghMinusGH=lh.amiValidated-ref.amiValidated;
        S.meanDeltaValidatedAMILoghMinusGH=mean(S.deltaValidatedAMILoghMinusGH,'omitnan');
        S.geometryL2=nan(nRep,1); S.geometryMaxAbs=nan(nRep,1);
        for r=1:nRep
            d=lhX{r}-refX{r};
            S.geometryL2(r)=norm(d,2);
            S.geometryMaxAbs(r)=max(abs(d));
        end
    else
        S.deltaValidatedAMILoghMinusGH=nan(nRep,1);
        S.meanDeltaValidatedAMILoghMinusGH=NaN;
        S.geometryL2=nan(nRep,1); S.geometryMaxAbs=nan(nRep,1);
    end
end


function out = stats_struct(raw)
    names=fieldnames(raw);
    out=struct();
    for k=1:numel(names)
        a=raw.(names{k});
        out.mean.(names{k})=mean(a,'omitnan');
        out.std.(names{k})=std(a,0,'omitnan');
        out.min.(names{k})=min(a,[],'omitnan');
        out.max.(names{k})=max(a,[],'omitnan');
    end
end


function figureFiles = make_plots(bundle,refCases,loghCases,figDir,saveFigures)
    E=bundle.evaluatorBenchmark; S=bundle.summary;
    figureFiles={}; figs={}; names={};

    f1=figure('Name','I.4 evaluator accuracy and speed','Color','w');
    tl=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('Same-constellation evaluator benchmark: GH%d vs log-h',bundle.methods.reference.order));
    nexttile;
    bar([abs(E.ghError(:)) abs(E.loghError(:))]); grid on;
    xlabel('Historical restart-winner candidate'); ylabel('|fast - validator| [bits/symbol]');
    legend(sprintf('GH%d',bundle.methods.reference.order),'log-h','Location','best');
    nexttile;
    bar([E.ghRuntimeMedian(:) E.loghRuntimeMedian(:)]); grid on;
    xlabel('Historical restart-winner candidate'); ylabel('Median evaluator runtime [s]');
    legend(sprintf('GH%d',bundle.methods.reference.order),'log-h','Location','best');
    figs{end+1}=f1; names{end+1}='evaluator_accuracy_speed';

    if ~isempty(loghCases) && all(~cellfun(@isempty,loghCases))
        f2=figure('Name','I.4 optimizer validated performance','Color','w');
        tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(tl,'End-to-end optimizer performance (validated)');
        nexttile;
        errorbar([1 2],[S.reference.mean.amiValidated S.logh.mean.amiValidated], ...
            [S.reference.std.amiValidated S.logh.std.amiValidated],'o','LineWidth',1.3);
        grid on; xlim([0.5 2.5]); set(gca,'XTick',[1 2],'XTickLabel',{'GH reference','log-h'});
        ylabel('Validated AMI [bits/symbol]');
        nexttile;
        errorbar([1 2],[S.reference.mean.gainBits S.logh.mean.gainBits], ...
            [S.reference.std.gainBits S.logh.std.gainBits],'o','LineWidth',1.3);
        grid on; xlim([0.5 2.5]); set(gca,'XTick',[1 2],'XTickLabel',{'GH reference','log-h'});
        ylabel('Shaping gain [bits/symbol]');
        figs{end+1}=f2; names{end+1}='optimizer_validated_performance';

        f3=figure('Name','I.4 convergence and geometry','Color','w');
        tl=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
        title(tl,'Convergence trajectories and validated winner geometry');
        nexttile; hold on;
        plot_mean_history(refCases,bundle.sa.historyEvery,'GH reference');
        plot_mean_history(loghCases,bundle.sa.historyEvery,'log-h');
        grid on; xlabel('SA iteration'); ylabel('Best-so-far fast AMI'); legend('Location','southeast');
        nexttile; hold on;
        M=bundle.case.M;
        for r=1:numel(refCases)
            plot(1:M,refCases{r}.winner.x(:),'o-','DisplayName',sprintf('GH rep%d',bundle.case.replicateVec(r)));
            plot(1:M,loghCases{r}.winner.x(:),'x--','DisplayName',sprintf('log-h rep%d',bundle.case.replicateVec(r)));
        end
        grid on; xlabel('Constellation index'); ylabel('Optical intensity'); legend('Location','best');
        figs{end+1}=f3; names{end+1}='convergence_geometry';

        f4=figure('Name','I.4 optimizer computational cost','Color','w');
        tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        title(tl,'Optimizer computational cost and convergence');
        methodLabels={'GH reference','log-h'};
        nexttile;
        bar([S.reference.mean.wallSeconds S.logh.mean.wallSeconds]); grid on;
        set(gca,'XTick',[1 2],'XTickLabel',methodLabels); ylabel('Wall time per replicate [s]');
        nexttile;
        bar([S.reference.mean.evalPerSec S.logh.mean.evalPerSec]); grid on;
        set(gca,'XTick',[1 2],'XTickLabel',methodLabels); ylabel('Fast objective eval/s');
        nexttile;
        bar([S.reference.mean.conv99 S.logh.mean.conv99]); grid on;
        set(gca,'XTick',[1 2],'XTickLabel',methodLabels); ylabel('Mean k_{99}');
        nexttile;
        bar([S.reference.mean.maxValidationGap S.logh.mean.maxValidationGap]); grid on;
        set(gca,'XTick',[1 2],'XTickLabel',methodLabels); ylabel('Mean max |fast-validator|');
        figs{end+1}=f4; names{end+1}='optimizer_cost_convergence';
    end

    if saveFigures
        if ~exist(figDir,'dir'), mkdir(figDir); end
        for k=1:numel(figs)
            png=fullfile(figDir,[names{k} '.png']);
            fig=fullfile(figDir,[names{k} '.fig']);
            exportgraphics(figs{k},png,'Resolution',220); savefig(figs{k},fig);
            figureFiles{end+1}=png; %#ok<AGROW>
            figureFiles{end+1}=fig; %#ok<AGROW>
        end
    end
end


function plot_mean_history(cases,historyEvery,label)
    maxIter=0;
    for r=1:numel(cases), maxIter=max(maxIter,cases{r}.winner.history.iter(end)); end
    common=unique([0:historyEvery:maxIter maxIter]);
    curves=nan(numel(cases),numel(common));
    for r=1:numel(cases)
        h=cases{r}.winner.history;
        curves(r,:)=interp1(h.iter,h.bestMI,common,'previous','extrap');
    end
    plot(common,mean(curves,1,'omitnan'),'LineWidth',1.4,'DisplayName',label);
end


function print_summary(bundle)
    S=bundle.summary;
    fprintf('\nCOMMIT I.4 SUMMARY\n');
    fprintf('Case: M=%d, SNR=%.2f dB, sigma_X^2=%.4f, d_min=%.4f\n', ...
        bundle.case.M,bundle.case.SNR_dB,bundle.case.sigma_X_sq,bundle.case.minGap);
    fprintf('Evaluator-only (%d historical restart winners):\n',numel(bundle.evaluatorBenchmark.validatorAMI));
    fprintf('  GH mean/max abs error   : %.3e / %.3e\n',S.evaluator.meanAbsErrorGH,S.evaluator.maxAbsErrorGH);
    fprintf('  log-h mean/max abs error: %.3e / %.3e\n',S.evaluator.meanAbsErrorLogh,S.evaluator.maxAbsErrorLogh);
    fprintf('  median runtime GH/log-h : %.4fs / %.4fs | speedup %.3fx\n', ...
        S.evaluator.medianRuntimeGH,S.evaluator.medianRuntimeLogh,S.evaluator.medianSpeedupGHOverLogh);
    fprintf('  max GH reproduction err : %.3e\n',S.evaluator.maxGHReproductionError);
    fprintf('  max validator repro err  : %.3e\n',S.evaluator.maxValidatorReproductionError);

    if all(isfinite(S.logh.raw.amiValidated))
        fprintf('\nOptimizer-level results:\n');
        fprintf(' rep | GH val AMI | log-h val AMI | delta(logh-GH) | GH k99 | LH k99 | GH eval/s | LH eval/s\n');
        fprintf('%s\n',repmat('-',1,100));
        for r=1:numel(S.replicateVec)
            fprintf('%4d | %.9f | %.9f | %+.3e | %6.0f | %6.0f | %9.3f | %9.3f\n', ...
                S.replicateVec(r),S.reference.raw.amiValidated(r),S.logh.raw.amiValidated(r), ...
                S.deltaValidatedAMILoghMinusGH(r),S.reference.raw.conv99(r),S.logh.raw.conv99(r), ...
                S.reference.raw.evalPerSec(r),S.logh.raw.evalPerSec(r));
        end
        fprintf('Mean validated AMI: GH=%.9f +/- %.2e | log-h=%.9f +/- %.2e\n', ...
            S.reference.mean.amiValidated,S.reference.std.amiValidated, ...
            S.logh.mean.amiValidated,S.logh.std.amiValidated);
        fprintf('Mean shaping gain : GH=%.9f | log-h=%.9f\n',S.reference.mean.gainBits,S.logh.mean.gainBits);
        fprintf('Mean delta log-h-GH validated AMI: %+.3e bits/symbol\n',S.meanDeltaValidatedAMILoghMinusGH);
        fprintf('Mean k99: GH=%.1f | log-h=%.1f\n',S.reference.mean.conv99,S.logh.mean.conv99);
        fprintf('Mean eval/s: GH=%.3f | log-h=%.3f\n',S.reference.mean.evalPerSec,S.logh.mean.evalPerSec);
        fprintf('Log-h fast/validated selection changes: %.0f of %d replicates\n', ...
            sum(S.logh.raw.selectionChanged,'omitnan'),numel(S.replicateVec));
        fprintf('Budget comparable to stored GH reference: %d\n',S.budgetComparable);
    end
end


function [src,path,dirPath] = load_reference_bundle(source)
    path=char(string(source));
    if ~isfile(path)
        error('sim_fast_evaluator_optimizer_benchmark:MissingSource','Source file not found: %s',path);
    end
    s=load(path);
    if ~isfield(s,'bundle') || ~isstruct(s.bundle)
        error('sim_fast_evaluator_optimizer_benchmark:BadSource','Source MAT must contain struct variable bundle.');
    end
    src=s.bundle;
    dirPath=fileparts(path);
end


function validate_reference_bundle(src)
    required={'experiment','case','settings','baselineSA','baselinePAM','summary','caseFiles'};
    for k=1:numel(required)
        if ~isfield(src,required{k})
            error('sim_fast_evaluator_optimizer_benchmark:BadSource','Source bundle missing field %s.',required{k});
        end
    end
    if ~strcmp(src.experiment,'sa_hyperparameter_sensitivity')
        error('sim_fast_evaluator_optimizer_benchmark:BadExperiment', ...
            'Expected Commit-G sa_hyperparameter_sensitivity source, got %s.',src.experiment);
    end
end


function [c,path] = load_reference_case(src,sourceDir,baselineIdx,repIdx)
    path='';
    try
        path=char(string(src.caseFiles{baselineIdx,repIdx}));
    catch
    end
    if isempty(path) || ~isfile(path)
        fname=sprintf('case_def%02d_rep%02d.mat',baselineIdx,repIdx);
        path=fullfile(sourceDir,'cases',fname);
    end
    if ~isfile(path)
        error('sim_fast_evaluator_optimizer_benchmark:MissingCase', ...
            'Could not locate reference baseline case for replicate %d.',repIdx);
    end
    s=load(path);
    if ~isfield(s,'c') || ~isstruct(s.c)
        error('sim_fast_evaluator_optimizer_benchmark:BadCase','Case file %s does not contain struct c.',path);
    end
    c=s.c;
end


function validate_reference_case(c,repIdx,M,SNR,sig,minGap)
    if ~isfield(c,'case') || ~isfield(c,'starts') || ~isfield(c,'winner')
        error('sim_fast_evaluator_optimizer_benchmark:BadCase','Reference case is missing required fields.');
    end
    if c.case.replicate~=repIdx || c.case.M~=M || ...
            abs(c.case.SNR_dB-SNR)>1e-12 || abs(c.case.sigma_X_sq-sig)>1e-12 || ...
            abs(c.case.minGap-minGap)>1e-12
        error('sim_fast_evaluator_optimizer_benchmark:CaseMismatch', ...
            'Reference case metadata does not match source bundle.');
    end
end


function id = source_run_id(src)
    id='';
    if isfield(src,'run') && isfield(src.run,'id'), id=src.run.id; end
end


function k = convergence_iteration(history,initialMI,bestMI,fraction)
    improvement=bestMI-initialMI;
    if ~(isfinite(improvement) && improvement>0), k=0; return; end
    target=initialMI+fraction*improvement;
    idx=find(history.bestMI>=target,1,'first');
    if isempty(idx), k=history.iter(end); else, k=history.iter(idx); end
end


function tf = positive_scalar(v)
    tf=isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v>0;
end
function tf = positive_integer(v)
    tf=positive_scalar(v) && mod(v,1)==0;
end
function tf = nonnegative_integer(v)
    tf=isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v>=0 && mod(v,1)==0;
end
function tf = positive_integer_vector(v)
    tf=isnumeric(v) && isvector(v) && ~isempty(v) && all(isfinite(v)) && all(v>=1) && all(mod(v,1)==0);
end


function token = num_token(v,nDec)
    token=sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token=strrep(token,'-','m'); token=strrep(token,'+','p'); token=strrep(token,'.','p');
end
function token = sanitize_run_label(label)
    token=char(string(label)); token=strtrim(token);
    if isempty(token), token='run'; end
    token=regexprep(token,'[^A-Za-z0-9_-]+','_'); token=regexprep(token,'_+','_');
    token=regexprep(token,'^_+|_+$',''); if isempty(token), token='run'; end
end
function [outDir,runId] = make_unique_run_directory(parentDir,requestedRunId)
    if ~exist(parentDir,'dir'), mkdir(parentDir); end
    runId=requestedRunId; outDir=fullfile(parentDir,runId); suffix=1;
    while exist(outDir,'dir')
        runId=sprintf('%s_%02d',requestedRunId,suffix); outDir=fullfile(parentDir,runId); suffix=suffix+1;
    end
    mkdir(outDir);
end
