function bundle = sim_sa_hyperparam_sensitivity(varargin)
%SIM_SA_HYPERPARAM_SENSITIVITY  SA hyperparameter sensitivity experiment.
%
% Commit G: reviewer-response experiment for SA validation/reproducibility.
% A representative FSO case is held fixed while the principal SA
% hyperparameters are varied one factor at a time (OFAT).
%
% Default representative case:
%   M = 32, SNR = 20 dB, sigma_X^2 = 0.1, d_min = 0.01
%   GH order = 400
%
% Default sweeps:
%   T0               = [0.1 0.2 0.4 0.8]
%   Tf               = [1e-4 1e-3 1e-2]
%   baseStd0         = [0.075 0.15 0.30]
%   itersPerTemp     = [25 50 100]
%
% Baseline SA configuration:
%   T0=0.4, Tf=1e-3, baseStd0=0.15, itersPerTemp=50.
%
% The baseline configuration is evaluated only once per replicate and is
% reused in all four sweep views.  For a fixed replicate, the same master
% seed is reused across SA settings.  Thus restart initializations and the
% initial random stream are paired, although trajectories may diverge after
% different accept/reject decisions.
%
% Reported metrics include independently validated AMI, shaping gain,
% convergence iterations k95/k99, runtime, fast-evaluator throughput,
% acceptance rate, restart variability, and fast-vs-validation discrepancy.
%
% Figures:
%   1) validated AMI and shaping gain versus each hyperparameter;
%   2) k99 and wall-clock runtime versus each hyperparameter;
%   3) best-so-far fast-AMI convergence trajectories;
%   4) acceptance rate and evaluator throughput.

    p = inputParser;
    p.FunctionName = mfilename;

    addParameter(p,'M',32,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=2 && mod(v,1)==0);
    addParameter(p,'SNR_dB',20,@(v) isnumeric(v) && isscalar(v) && isfinite(v));
    addParameter(p,'sigma_X_sq',0.1,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p,'P_avg',1,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p,'minGap',0.01,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p,'ghN_h',400,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=2 && mod(v,1)==0);

    addParameter(p,'saMaxIter',2000,@positive_integer);
    addParameter(p,'saNStarts',4,@(v) positive_integer(v) && v>=2);
    addParameter(p,'nReplicates',2,@positive_integer);
    addParameter(p,'baseSeed',20260928,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p,'historyEvery',50,@positive_integer);

    addParameter(p,'BaselineT0',0.4,@positive_scalar);
    addParameter(p,'BaselineTf',1e-3,@positive_scalar);
    addParameter(p,'BaselineBaseStd0',0.15,@positive_scalar);
    addParameter(p,'BaselineItersPerTemp',50,@positive_integer);

    addParameter(p,'T0Vec',[0.1 0.2 0.4 0.8],@positive_vector);
    addParameter(p,'TfVec',[1e-4 1e-3 1e-2],@positive_vector);
    addParameter(p,'BaseStd0Vec',[0.075 0.15 0.30],@positive_vector);
    addParameter(p,'ItersPerTempVec',[25 50 100],@positive_integer_vector);

    addParameter(p,'UseParallelCases',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'MakePlots',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'SaveFigures',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v) || isstring(v));
    parse(p,varargin{:});
    o = p.Results;

    validate_baseline_membership(o);

    requiredMean = (double(o.M)-1)*double(o.minGap)/2;
    if requiredMean > double(o.P_avg) + 100*eps(max(1,double(o.P_avg)))
        error('sim_sa_hyperparam_sensitivity:InfeasibleGap', ...
            'M=%d and d_min=%.6g require mean power %.6g > P_avg=%.6g.', ...
            o.M,o.minGap,requiredMean,o.P_avg);
    end

    base = struct( ...
        'id',1, ...
        'name','baseline', ...
        'T0',double(o.BaselineT0), ...
        'Tf',double(o.BaselineTf), ...
        'baseStd0',double(o.BaselineBaseStd0), ...
        'itersPerTemp',double(o.BaselineItersPerTemp));

    defs = base;

    T0Vec  = unique(double(o.T0Vec(:).'),'stable');
    TfVec  = unique(double(o.TfVec(:).'),'stable');
    StdVec = unique(double(o.BaseStd0Vec(:).'),'stable');
    BlkVec = unique(double(o.ItersPerTempVec(:).'),'stable');

    for v = T0Vec
        d = strip_identity(base); d.T0 = v;
        defs = add_unique_def(defs,d,sprintf('T0=%.6g',v));
    end
    for v = TfVec
        d = strip_identity(base); d.Tf = v;
        defs = add_unique_def(defs,d,sprintf('Tf=%.6g',v));
    end
    for v = StdVec
        d = strip_identity(base); d.baseStd0 = v;
        defs = add_unique_def(defs,d,sprintf('baseStd0=%.6g',v));
    end
    for v = BlkVec
        d = strip_identity(base); d.itersPerTemp = v;
        defs = add_unique_def(defs,d,sprintf('itersPerTemp=%d',v));
    end

    nDefs = numel(defs);
    nRep  = double(o.nReplicates);

    % IMPORTANT: initialize from the first real structure.  MATLAB does not
    % permit assigning a field-bearing structure into struct([]) on all
    % releases (the original Commit-G smoke test exposed this on Windows).
    sweeps = make_sweep('T0','Initial temperature T_0',T0Vec,defs,base);
    sweeps(2) = make_sweep('Tf','Final temperature T_f',TfVec,defs,base);
    sweeps(3) = make_sweep('baseStd0','Initial proposal std',StdVec,defs,base);
    sweeps(4) = make_sweep('itersPerTemp','Iterations per temperature block',BlkVec,defs,base);

    resultsRoot = char(o.ResultsRoot);
    conditionTag = sprintf('M%d_SNR%s_sig%s_gap%s', ...
        o.M,num_token(o.SNR_dB,2),num_token(o.sigma_X_sq,4),num_token(o.minGap,4));
    conditionDir = fullfile(resultsRoot,'sa_hyperparameter_sensitivity',conditionTag);

    runMeta = fso_result_utils.run_metadata(mfilename);
    labelToken = sanitize_run_label(o.RunLabel);
    configToken = sprintf('GH%d_I%d_S%d_R%d',o.ghN_h,o.saMaxIter,o.saNStarts,nRep);
    utcToken = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requestedRunId = sprintf('%s_%s_%s',labelToken,configToken,utcToken);
    [outDir,runId] = make_unique_run_directory(conditionDir,requestedRunId);

    caseDir = fullfile(outDir,'cases');
    figDir  = fullfile(outDir,'figures');
    mkdir(caseDir);
    if o.MakePlots && o.SaveFigures, mkdir(figDir); end

    runMeta.runId = runId;
    runMeta.runLabel = char(string(o.RunLabel));
    runMeta.outputDirectory = outDir;
    runMeta.conditionTag = conditionTag;

    xMaxBound = fso_result_utils.feasible_xmax_bound(o.M,o.P_avg,o.minGap,true);
    cfgBase = build_fso_config(o.M,o.P_avg,o.SNR_dB,o.sigma_X_sq, ...
        'ghN_h',o.ghN_h,'xMaxBound',xMaxBound, ...
        'saMaxIter',o.saMaxIter,'saNStarts',o.saNStarts, ...
        'saUseParallel',false,'minGap',o.minGap,'pinZero',true, ...
        'seedInit',1,'logEvery',0,'historyEvery',o.historyEvery);

    xPAM = define_constellation(o.M,o.P_avg,true,"mean");
    AMI_functions.assert_constellation_feasible(xPAM,cfgBase,1e-10);
    amiPAMFast = cfgBase.AMI_Evaluator(xPAM);
    amiPAMValidated = cfgBase.AMI_Validator(xPAM);

    fprintf('\n============================================================\n');
    fprintf('IEEE revision SA hyperparameter sensitivity experiment\n');
    fprintf('Run ID: %s\n',runId);
    if strlength(string(o.RunLabel)) > 0
        fprintf('Run label: %s\n',char(string(o.RunLabel)));
    end
    fprintf('Representative case: M=%d | SNR=%.2f dB | sigma_X^2=%.4f | d_min=%.4f\n', ...
        o.M,o.SNR_dB,o.sigma_X_sq,o.minGap);
    fprintf('GH=%d | SA=%d restarts x %d iterations | replicates=%d\n', ...
        o.ghN_h,o.saNStarts,o.saMaxIter,nRep);
    fprintf('Baseline SA: T0=%.4g | Tf=%.4g | step0=%.4g | block=%d\n', ...
        base.T0,base.Tf,base.baseStd0,base.itersPerTemp);
    fprintf('T0 sweep: %s\n',mat2str(T0Vec));
    fprintf('Tf sweep: %s\n',mat2str(TfVec));
    fprintf('step0 sweep: %s\n',mat2str(StdVec));
    fprintf('block sweep: %s\n',mat2str(BlkVec));
    fprintf('Unique SA settings=%d | total cases=%d\n',nDefs,nDefs*nRep);
    fprintf('Paired master seeds across SA settings: enabled\n');
    fprintf('Uniform-PAM validated AMI=%.9f\n',amiPAMValidated);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    [iDGrid,iRGrid] = ndgrid(1:nDefs,1:nRep);
    taskD = iDGrid(:);
    taskR = iRGrid(:);
    nCases = numel(taskD);

    caseResults = cell(nCases,1);
    caseFiles   = cell(nCases,1);

    doPar = logical(o.UseParallelCases);
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
        fprintf('Outer case parallelism: %d process workers\n',pool.NumWorkers);

        Q = parallel.pool.DataQueue;
        doneCount = 0;
        afterEach(Q,@(~) progress_tick());

        parfor k = 1:nCases
            [caseResults{k},caseFiles{k}] = run_sa_case( ...
                defs(taskD(k)),taskR(k),o,runMeta,caseDir,xPAM, ...
                amiPAMFast,amiPAMValidated,xMaxBound);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k = 1:nCases
            [caseResults{k},caseFiles{k}] = run_sa_case( ...
                defs(taskD(k)),taskR(k),o,runMeta,caseDir,xPAM, ...
                amiPAMFast,amiPAMValidated,xMaxBound);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    caseGrid = reshape(caseResults,nDefs,nRep);
    fileGrid = reshape(caseFiles,nDefs,nRep);
    summary = build_summary(caseGrid,defs,nRep,amiPAMValidated,o.saMaxIter);

    bundle = struct();
    bundle.meta = runMeta;
    bundle.experiment = 'sa_hyperparameter_sensitivity';
    bundle.run = struct('id',runId,'label',char(string(o.RunLabel)), ...
        'conditionTag',conditionTag,'outputDirectory',outDir);
    bundle.case = struct('M',double(o.M),'P_avg',double(o.P_avg), ...
        'SNR_dB',double(o.SNR_dB),'sigma_X_sq',double(o.sigma_X_sq), ...
        'minGap',double(o.minGap),'ghN_h',double(o.ghN_h));
    bundle.baselinePAM = struct('x',xPAM(:),'amiFast',amiPAMFast, ...
        'amiValidated',amiPAMValidated);
    bundle.baselineSA = base;
    bundle.settings = struct('saMaxIter',double(o.saMaxIter), ...
        'saNStarts',double(o.saNStarts),'nReplicates',nRep, ...
        'baseSeed',double(o.baseSeed),'historyEvery',double(o.historyEvery), ...
        'pairedMasterSeedsAcrossSettings',true,'outerParallelism',doPar);
    bundle.definitions = defs;
    bundle.sweeps = sweeps;
    bundle.caseFiles = fileGrid;
    bundle.summary = summary;
    bundle.figureFiles = {};

    if o.MakePlots
        bundle.figureFiles = make_plots(bundle,caseGrid,figDir,logical(o.SaveFigures));
    end

    bundlePath = fullfile(outDir,'saHyperparameterSensitivity.mat');
    save(bundlePath,'bundle','-v7.3');

    fprintf('\nSaved summary bundle: %s\n',bundlePath);
    print_summary(bundle);

    function progress_tick()
        doneCount = doneCount + 1;
        fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function [c,path] = run_sa_case(def,repIdx,o,runMeta,caseDir,xPAM,amiPAMFast,amiPAMValidated,xMaxBound)
    replicateBaseSeed = double(o.baseSeed) + (repIdx-1)*10000019;
    seed = fso_result_utils.case_seed(replicateBaseSeed,o.M,o.SNR_dB,o.sigma_X_sq);

    cfg = build_fso_config(o.M,o.P_avg,o.SNR_dB,o.sigma_X_sq, ...
        'ghN_h',o.ghN_h,'xMaxBound',xMaxBound, ...
        'saMaxIter',o.saMaxIter,'saNStarts',o.saNStarts, ...
        'saUseParallel',false,'minGap',o.minGap,'pinZero',true, ...
        'seedInit',seed,'logEvery',0,'historyEvery',o.historyEvery);

    cfg.SA.T0 = def.T0;
    cfg.SA.Tf = def.Tf;
    cfg.SA.baseStd0 = def.baseStd0;
    cfg.SA.itersPerTemp = def.itersPerTemp;
    cfg.SA.nBlocks = ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate = exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);

    if ~(cfg.SA.T0 > cfg.SA.Tf)
        error('sim_sa_hyperparam_sensitivity:TemperatureOrder', ...
            'Configuration %s has T0 <= Tf.',def.name);
    end

    AMI_functions.assert_constellation_feasible(xPAM,cfg,1e-10);

    label = sprintf('[SAparam|%s|rep%d]',def.name,repIdx);
    tCase = tic;
    [saOut,starts] = sa_multistart(cfg,xPAM,label);
    wallRuntime = toc(tCase);

    win = starts(saOut.bestStartValidated);
    k95 = convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.95);
    k99 = convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.99);

    vals = [starts.bestMIValidated];
    totalFastEvaluations = double(o.saNStarts)*double(o.saMaxIter);
    evalThroughput = totalFastEvaluations/max(saOut.totalSARuntime,eps);

    c = struct();
    c.meta = runMeta;
    c.case = struct('definitionId',def.id,'definitionName',def.name, ...
        'replicate',repIdx,'seed',seed,'M',double(o.M), ...
        'SNR_dB',double(o.SNR_dB),'sigma_X_sq',double(o.sigma_X_sq), ...
        'minGap',double(o.minGap));
    c.sa = struct('T0',def.T0,'Tf',def.Tf,'baseStd0',def.baseStd0, ...
        'itersPerTemp',def.itersPerTemp,'coolingRate',cfg.SA.coolingRate, ...
        'maxIter',cfg.SA.maxIter,'nStarts',cfg.SA.nStarts, ...
        'targetAccLo',cfg.SA.targetAccLo,'targetAccHi',cfg.SA.targetAccHi, ...
        'baseStdMin',cfg.SA.baseStdMin,'baseStdMax',cfg.SA.baseStdMax, ...
        'baseStdGrow',cfg.SA.baseStdGrow,'baseStdShrink',cfg.SA.baseStdShrink);
    c.config = fso_result_utils.config_snapshot(cfg);
    c.baseline = struct('x',xPAM(:),'amiFast',amiPAMFast, ...
        'amiValidated',amiPAMValidated);
    c.starts = starts;
    c.optimizer = saOut;
    if isfield(c.optimizer,'bestRun'), c.optimizer = rmfield(c.optimizer,'bestRun'); end
    c.winner = struct('startIndex',saOut.bestStartValidated, ...
        'x',saOut.bestXValidated(:), ...
        'amiFast',saOut.validatedWinnerFastMI, ...
        'amiValidated',saOut.bestMIValidated, ...
        'gainBits',saOut.bestMIValidated-amiPAMValidated, ...
        'acceptanceRate',win.acceptanceRate, ...
        'convergenceIter95',k95,'convergenceIter99',k99, ...
        'history',win.history);
    c.repeatability = struct('restartValidatedMean',mean(vals,'omitnan'), ...
        'restartValidatedStd',std(vals,0,'omitnan'), ...
        'restartValidatedRange',max(vals)-min(vals));
    c.validation = struct('selectionChanged',saOut.selectionChanged, ...
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap, ...
        'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    c.runtime = struct('wallSeconds',wallRuntime, ...
        'saSecondsSummed',saOut.totalSARuntime, ...
        'validationSeconds',saOut.validationTotalRuntime, ...
        'fastEvaluations',totalFastEvaluations, ...
        'fastEvaluationsPerSecond',evalThroughput);

    fname = sprintf('case_def%02d_rep%02d.mat',def.id,repIdx);
    path = fullfile(caseDir,fname);
    save(path,'c','-v7.3');
end


function summary = build_summary(caseGrid,defs,nRep,amiPAMValidated,maxIter)
    nDefs = numel(defs);
    metricNames = {'amiValidated','gainBits','conv95','conv99','wallSeconds', ...
        'evalPerSec','acceptanceRate','restartStd','maxValidationGap','selectionChanged'};

    raw = struct();
    for f = 1:numel(metricNames)
        raw.(metricNames{f}) = nan(nDefs,nRep);
    end

    for iD = 1:nDefs
        for iR = 1:nRep
            c = caseGrid{iD,iR};
            raw.amiValidated(iD,iR) = c.winner.amiValidated;
            raw.gainBits(iD,iR) = c.winner.gainBits;
            raw.conv95(iD,iR) = c.winner.convergenceIter95;
            raw.conv99(iD,iR) = c.winner.convergenceIter99;
            raw.wallSeconds(iD,iR) = c.runtime.wallSeconds;
            raw.evalPerSec(iD,iR) = c.runtime.fastEvaluationsPerSecond;
            raw.acceptanceRate(iD,iR) = c.winner.acceptanceRate;
            raw.restartStd(iD,iR) = c.repeatability.restartValidatedStd;
            raw.maxValidationGap(iD,iR) = c.validation.maxAbsFastValidationGap;
            raw.selectionChanged(iD,iR) = double(c.validation.selectionChanged);
        end
    end

    summary = struct();
    summary.definitionNames = {defs.name};
    summary.raw = raw;
    summary.baselinePAMValidated = amiPAMValidated;
    summary.maxIter = maxIter;

    names = fieldnames(raw);
    for k = 1:numel(names)
        a = raw.(names{k});
        summary.mean.(names{k}) = mean(a,2,'omitnan');
        summary.std.(names{k})  = std(a,0,2,'omitnan');
        summary.min.(names{k})  = min(a,[],2,'omitnan');
        summary.max.(names{k})  = max(a,[],2,'omitnan');
    end

    baselineIdx = find(strcmp({defs.name},'baseline'),1);
    summary.baselineDefinitionIndex = baselineIdx;
    summary.deltaAMIvsBaseline = raw.amiValidated - raw.amiValidated(baselineIdx,:);
    summary.meanDeltaAMIvsBaseline = mean(summary.deltaAMIvsBaseline,2,'omitnan');
end


function figureFiles = make_plots(bundle,caseGrid,figDir,saveFigures)
    sweeps = bundle.sweeps;
    S = bundle.summary;
    nRep = bundle.settings.nReplicates;
    figureFiles = {};

    f1 = figure('Name','SA hyperparameters - validated performance','Color','w');
    tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('SA sensitivity: M=%d, SNR=%.1f dB, \\sigma_X^2=%.3f, d_{min}=%.3f', ...
        bundle.case.M,bundle.case.SNR_dB,bundle.case.sigma_X_sq,bundle.case.minGap));
    for s = 1:numel(sweeps)
        idx = sweeps(s).definitionIndices;
        x = sweeps(s).values;
        nexttile;
        errorbar(x,S.mean.amiValidated(idx),S.std.amiValidated(idx),'o-','LineWidth',1.3);
        grid on; xlabel(sweeps(s).axisLabel); ylabel('Validated AMI [bits/symbol]');
        apply_log_axis_if_needed(sweeps(s).field);
        title(sprintf('%s: AMI',sweeps(s).displayName));
        nexttile;
        errorbar(x,S.mean.gainBits(idx),S.std.gainBits(idx),'o-','LineWidth',1.3);
        grid on; xlabel(sweeps(s).axisLabel); ylabel('Shaping gain [bits/symbol]');
        apply_log_axis_if_needed(sweeps(s).field);
        title(sprintf('%s: gain',sweeps(s).displayName));
    end

    f2 = figure('Name','SA hyperparameters - convergence and runtime','Color','w');
    tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
    title(tl,'Convergence speed and computational cost');
    for s = 1:numel(sweeps)
        idx = sweeps(s).definitionIndices;
        x = sweeps(s).values;
        nexttile;
        errorbar(x,S.mean.conv99(idx),S.std.conv99(idx),'o-','LineWidth',1.3);
        hold on; yline(bundle.settings.saMaxIter,'--','Iteration budget'); hold off;
        grid on; xlabel(sweeps(s).axisLabel); ylabel('Iteration reaching 99% improvement');
        apply_log_axis_if_needed(sweeps(s).field);
        title(sprintf('%s: k_{99}',sweeps(s).displayName));
        nexttile;
        errorbar(x,S.mean.wallSeconds(idx),S.std.wallSeconds(idx),'o-','LineWidth',1.3);
        grid on; xlabel(sweeps(s).axisLabel); ylabel('Wall time per case [s]');
        apply_log_axis_if_needed(sweeps(s).field);
        title(sprintf('%s: runtime',sweeps(s).displayName));
    end

    f3 = figure('Name','SA hyperparameters - convergence trajectories','Color','w');
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(tl,'Best-so-far fast-objective convergence trajectories');
    commonIter = unique([0:bundle.settings.historyEvery:bundle.settings.saMaxIter, ...
        bundle.settings.saMaxIter]);
    for s = 1:numel(sweeps)
        nexttile; hold on;
        idx = sweeps(s).definitionIndices;
        for j = 1:numel(idx)
            curves = nan(nRep,numel(commonIter));
            for r = 1:nRep
                h = caseGrid{idx(j),r}.winner.history;
                curves(r,:) = interp1(h.iter,h.bestMI,commonIter,'previous','extrap');
            end
            plot(commonIter,mean(curves,1,'omitnan'),'LineWidth',1.3, ...
                'DisplayName',sprintf('%s=%g',sweeps(s).field,sweeps(s).values(j)));
        end
        grid on; xlabel('SA iteration'); ylabel('Best-so-far fast AMI');
        title(sweeps(s).displayName); legend('Location','southeast');
    end

    f4 = figure('Name','SA hyperparameters - acceptance and throughput','Color','w');
    tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
    title(tl,'Search dynamics and evaluator throughput');
    for s = 1:numel(sweeps)
        idx = sweeps(s).definitionIndices;
        x = sweeps(s).values;
        nexttile;
        errorbar(x,S.mean.acceptanceRate(idx),S.std.acceptanceRate(idx),'o-','LineWidth',1.3);
        hold on; yline(0.20,'--'); yline(0.60,'--'); hold off;
        ylim([0 1]); grid on; xlabel(sweeps(s).axisLabel); ylabel('Winner acceptance rate');
        apply_log_axis_if_needed(sweeps(s).field);
        title(sprintf('%s: acceptance',sweeps(s).displayName));
        nexttile;
        errorbar(x,S.mean.evalPerSec(idx),S.std.evalPerSec(idx),'o-','LineWidth',1.3);
        grid on; xlabel(sweeps(s).axisLabel); ylabel('Fast objective eval/s');
        apply_log_axis_if_needed(sweeps(s).field);
        title(sprintf('%s: throughput',sweeps(s).displayName));
    end

    if saveFigures
        if ~exist(figDir,'dir'), mkdir(figDir); end
        figs = {f1,f2,f3,f4};
        names = {'validated_performance','convergence_runtime', ...
            'convergence_trajectories','acceptance_throughput'};
        for k = 1:numel(figs)
            pngPath = fullfile(figDir,[names{k} '.png']);
            figPath = fullfile(figDir,[names{k} '.fig']);
            exportgraphics(figs{k},pngPath,'Resolution',220);
            savefig(figs{k},figPath);
            figureFiles{end+1} = pngPath; %#ok<AGROW>
            figureFiles{end+1} = figPath; %#ok<AGROW>
        end
    end
end


function print_summary(bundle)
    fprintf('\nSA HYPERPARAMETER SENSITIVITY SUMMARY\n');
    fprintf('Representative case: M=%d, SNR=%.2f dB, sigma_X^2=%.4f, d_min=%.4f\n', ...
        bundle.case.M,bundle.case.SNR_dB,bundle.case.sigma_X_sq,bundle.case.minGap);
    fprintf('Validated Uniform-PAM AMI: %.9f\n',bundle.baselinePAM.amiValidated);

    S = bundle.summary;
    defs = bundle.definitions;
    for s = 1:numel(bundle.sweeps)
        sw = bundle.sweeps(s);
        fprintf('\n--- %s ---\n',sw.displayName);
        fprintf('%13s | validated AMI      | gain               | delta vs base | k99    | runtime[s] | eval/s   | accept | restart std | max F-V gap\n',sw.field);
        fprintf('%s\n',repmat('-',1,142));
        for j = 1:numel(sw.values)
            i = sw.definitionIndices(j);
            fprintf('%13.6g | %.6f +/- %.2e | %+.6f +/- %.2e | %+.3e | %6.1f | %10.1f | %8.2f | %.3f | %.3e | %.3e\n', ...
                sw.values(j),S.mean.amiValidated(i),S.std.amiValidated(i), ...
                S.mean.gainBits(i),S.std.gainBits(i),S.meanDeltaAMIvsBaseline(i), ...
                S.mean.conv99(i),S.mean.wallSeconds(i),S.mean.evalPerSec(i), ...
                S.mean.acceptanceRate(i),S.mean.restartStd(i),S.mean.maxValidationGap(i));
        end
    end

    [bestAMI,bestIdx] = max(S.mean.amiValidated);
    [fastConv,fastIdx] = min(S.mean.conv99);
    fprintf('\nBest observed mean validated AMI: %.9f (%s)\n',bestAMI,defs(bestIdx).name);
    fprintf('Fastest observed mean k99: %.1f iterations (%s)\n',fastConv,defs(fastIdx).name);
    fprintf('NOTE: descriptive sensitivity results do not prove globally optimal SA hyperparameters.\n');
end


function k = convergence_iteration(history,initialMI,bestMI,fraction)
    improvement = bestMI-initialMI;
    if ~(isfinite(improvement) && improvement>0)
        k = 0;
        return;
    end
    target = initialMI + fraction*improvement;
    idx = find(history.bestMI >= target,1,'first');
    if isempty(idx), k = history.iter(end); else, k = history.iter(idx); end
end


function sw = make_sweep(field,displayName,values,defs,base)
    idx = nan(size(values));
    for k = 1:numel(values)
        d = strip_identity(base);
        d.(field) = values(k);
        idx(k) = find_definition(defs,d,true);
    end
    sw = struct('field',field,'displayName',displayName,'values',values, ...
        'definitionIndices',idx,'axisLabel',axis_label(field));
end


function defs = add_unique_def(defs,d,name)
    idx = find_definition(defs,d,false);
    if ~isnan(idx), return; end
    entry = struct('id',numel(defs)+1,'name',name, ...
        'T0',double(d.T0),'Tf',double(d.Tf), ...
        'baseStd0',double(d.baseStd0),'itersPerTemp',double(d.itersPerTemp));
    defs(end+1) = entry;
end


function d = strip_identity(s)
    d = struct('T0',double(s.T0),'Tf',double(s.Tf), ...
        'baseStd0',double(s.baseStd0),'itersPerTemp',double(s.itersPerTemp));
end


function idx = find_definition(defs,d,mustExist)
    if nargin < 3, mustExist = true; end
    idx = NaN;
    for k = 1:numel(defs)
        if same_number(defs(k).T0,d.T0) && same_number(defs(k).Tf,d.Tf) && ...
                same_number(defs(k).baseStd0,d.baseStd0) && ...
                same_number(defs(k).itersPerTemp,d.itersPerTemp)
            idx = k;
            return;
        end
    end
    if mustExist
        error('sim_sa_hyperparam_sensitivity:MissingDefinition', ...
            'Internal sweep definition was not found.');
    end
end


function tf = same_number(a,b)
    tf = abs(double(a)-double(b)) <= 100*eps(max([1 abs(double(a)) abs(double(b))]));
end


function label = axis_label(field)
    switch field
        case 'T0', label = 'T_0';
        case 'Tf', label = 'T_f';
        case 'baseStd0', label = 'Initial proposal std';
        case 'itersPerTemp', label = 'Iterations per temperature block';
        otherwise, label = field;
    end
end


function apply_log_axis_if_needed(field)
    if strcmp(field,'T0') || strcmp(field,'Tf')
        set(gca,'XScale','log');
    end
end


function validate_baseline_membership(o)
    if ~any_close(o.T0Vec,o.BaselineT0)
        error('sim_sa_hyperparam_sensitivity:BaselineMissing','T0Vec must contain BaselineT0.');
    end
    if ~any_close(o.TfVec,o.BaselineTf)
        error('sim_sa_hyperparam_sensitivity:BaselineMissing','TfVec must contain BaselineTf.');
    end
    if ~any_close(o.BaseStd0Vec,o.BaselineBaseStd0)
        error('sim_sa_hyperparam_sensitivity:BaselineMissing','BaseStd0Vec must contain BaselineBaseStd0.');
    end
    if ~any_close(o.ItersPerTempVec,o.BaselineItersPerTemp)
        error('sim_sa_hyperparam_sensitivity:BaselineMissing','ItersPerTempVec must contain BaselineItersPerTemp.');
    end

    if any(double(o.TfVec) >= double(o.BaselineT0)) || ...
            any(double(o.T0Vec) <= double(o.BaselineTf))
        error('sim_sa_hyperparam_sensitivity:TemperatureOrder', ...
            'Every OFAT temperature setting must satisfy T0 > Tf.');
    end
end


function tf = any_close(vec,val)
    vec = double(vec(:)); val = double(val);
    tf = any(abs(vec-val) <= 100*eps(max(1,abs(val))));
end


function tf = positive_scalar(v)
    tf = isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v>0;
end

function tf = positive_integer(v)
    tf = positive_scalar(v) && mod(v,1)==0;
end

function tf = positive_vector(v)
    tf = isnumeric(v) && isvector(v) && ~isempty(v) && ...
        all(isfinite(v)) && all(isreal(v)) && all(v>0);
end

function tf = positive_integer_vector(v)
    tf = positive_vector(v) && all(mod(v,1)==0);
end


function token = num_token(v,nDec)
    token = sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token = strrep(token,'-','m');
    token = strrep(token,'+','p');
    token = strrep(token,'.','p');
end


function token = sanitize_run_label(label)
    token = char(string(label));
    token = strtrim(token);
    if isempty(token), token = 'run'; end
    token = regexprep(token,'[^A-Za-z0-9_-]+','_');
    token = regexprep(token,'_+','_');
    token = regexprep(token,'^_+|_+$','');
    if isempty(token), token = 'run'; end
end


function [outDir,runId] = make_unique_run_directory(parentDir,requestedRunId)
    if ~exist(parentDir,'dir'), mkdir(parentDir); end
    runId = requestedRunId;
    outDir = fullfile(parentDir,runId);
    suffix = 1;
    while exist(outDir,'dir')
        runId = sprintf('%s_%02d',requestedRunId,suffix);
        outDir = fullfile(parentDir,runId);
        suffix = suffix+1;
    end
    mkdir(outDir);
end
