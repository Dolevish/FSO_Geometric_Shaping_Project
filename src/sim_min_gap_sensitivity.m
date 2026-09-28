function bundle = sim_min_gap_sensitivity(varargin)
%SIM_MIN_GAP_SENSITIVITY  Reviewer-response experiment for d_min sensitivity.
%
%   bundle = sim_min_gap_sensitivity()
%   bundle = sim_min_gap_sensitivity('Name',Value,...)
%
% This experiment measures how the numerical minimum-spacing constraint
% affects the validated AMI, shaping gain, optimized geometry and SA
% robustness. It is designed to address the editor/Reviewer 1/Reviewer 4
% concern that d_min=0.05 may restrict the GS solution, especially for
% M=16 and M=32.
%
% Default research grid:
%   M                = [8 16 32]
%   SNR              = 20 dB
%   sigma_X^2        = 0.1
%   d_min            = [0 0.005 0.01 0.02 0.05]
%
% Name-value options:
%   'MVec'             constellation sizes (default [8 16 32])
%   'SNR_dB'           fixed SNR in dB (default 20)
%   'sigma_X_sq'       fixed normalized intensity variance (default 0.1)
%   'P_avg'            average optical intensity (default 1)
%   'MinGapVec'        imposed minimum gaps (default [0 .005 .01 .02 .05])
%   'ghN_h'            Gauss-Hermite order (default 40)
%   'saMaxIter'        iterations per restart (default 10000)
%   'saNStarts'        restarts per case (default 8)
%   'nReplicates'      independent multistart batches per gap (default 1)
%   'baseSeed'         deterministic experiment seed (default 20260928)
%   'historyEvery'     SA history cadence ([] -> temperature block)
%   'pinZero'          force first level to zero (default true)
%   'UseParallelCases' parallelize M/gap/replicate cases (default true)
%   'ComputeBER'       auxiliary BER diagnostic (default false)
%   'MakePlots'        create summary figures (default true)
%   'SaveFigures'      save PNG/FIG files when MakePlots=true (default true)
%   'RunLabel'         optional human-readable run label (default '')
%   'ResultsRoot'      output root (default repository/results/ieee_revision)
%
% RUN IDENTITY / NO-OVERWRITE GUARANTEE (Commit F.1):
% Every invocation is written to its own run directory. The directory name
% contains a sanitized RunLabel, numerical settings (GH/iterations/starts/
% replicates), and a UTC timestamp. If an unlikely name collision occurs, a
% numeric suffix is added. Existing experiment files are never overwritten.
%
% Methodological detail: for a fixed M and replicate, the SAME master seed
% is deliberately reused across all d_min values. Thus each d_min case uses
% paired stochastic restarts; only the feasible-set projection changes.
% Different replicates receive independent master seeds.
%
% Full per-case results include all restart histories from sa_multistart().
% The lightweight bundle stores summary arrays and file references.

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p,'MVec',[8 16 32],@(v) isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>=2) && all(mod(v,1)==0));
    addParameter(p,'SNR_dB',20,@(v) isnumeric(v) && isscalar(v) && isfinite(v));
    addParameter(p,'sigma_X_sq',0.1,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p,'P_avg',1,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p,'MinGapVec',[0 0.005 0.01 0.02 0.05],@(v) isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>=0));
    addParameter(p,'ghN_h',40,@(v) isnumeric(v) && isscalar(v) && v>=2 && mod(v,1)==0);
    addParameter(p,'saMaxIter',10000,@(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p,'saNStarts',8,@(v) isnumeric(v) && isscalar(v) && v>=2 && mod(v,1)==0);
    addParameter(p,'nReplicates',1,@(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p,'baseSeed',20260928,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p,'historyEvery',[],@(v) isempty(v) || (isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0));
    addParameter(p,'pinZero',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'UseParallelCases',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'ComputeBER',false,@(v) islogical(v) && isscalar(v));
    addParameter(p,'MakePlots',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'SaveFigures',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v) || isstring(v));
    parse(p,varargin{:});
    o = p.Results;

    MVec   = unique(double(o.MVec(:).'),'stable');
    gapVec = unique(double(o.MinGapVec(:).'),'stable');
    nM     = numel(MVec);
    nGap   = numel(gapVec);
    nRep   = double(o.nReplicates);

    % Fail before a long run if any requested d_min is infeasible.
    for iM = 1:nM
        requiredMean = (MVec(iM)-1).*gapVec/2;
        if any(requiredMean > o.P_avg + 100*eps(max(1,o.P_avg)))
            bad = find(requiredMean > o.P_avg + 100*eps(max(1,o.P_avg)),1,'first');
            error('sim_min_gap_sensitivity:InfeasibleGap', ...
                'M=%d and d_min=%.6g require mean power %.6g > P_avg=%.6g.', ...
                MVec(iM),gapVec(bad),requiredMean(bad),o.P_avg);
        end
    end

    % ------------------------------------------------------------------
    % Commit F.1: immutable run identity and isolated output directory.
    % ------------------------------------------------------------------
    resultsRoot = char(o.ResultsRoot);
    conditionTag = sprintf('SNR%s_sig%s',num_token(o.SNR_dB,2),num_token(o.sigma_X_sq,4));
    conditionDir = fullfile(resultsRoot,'min_gap_sensitivity',conditionTag);

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

    [iMGrid,iGGrid,iRGrid] = ndgrid(1:nM,1:nGap,1:nRep);
    taskM   = iMGrid(:);
    taskG   = iGGrid(:);
    taskR   = iRGrid(:);
    nCases  = numel(taskM);

    fprintf('\n============================================================\n');
    fprintf('IEEE revision minimum-spacing sensitivity experiment\n');
    fprintf('Run ID: %s\n',runId);
    if strlength(string(o.RunLabel)) > 0
        fprintf('Run label: %s\n',char(string(o.RunLabel)));
    end
    fprintf('M = %s | SNR=%.2f dB | sigma_X^2=%.4f | Pavg=%.4g\n', ...
        mat2str(MVec),o.SNR_dB,o.sigma_X_sq,o.P_avg);
    fprintf('d_min = %s\n',mat2str(gapVec));
    fprintf('Replicates=%d | SA=%d restarts x %d iterations | GH=%d\n', ...
        nRep,o.saNStarts,o.saMaxIter,o.ghN_h);
    fprintf('Paired seeds across d_min: enabled\n');
    fprintf('Total optimization cases: %d\n',nCases);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

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
            [caseResults{k},caseFiles{k}] = run_gap_case( ...
                MVec(taskM(k)),gapVec(taskG(k)),taskR(k),o,runMeta,caseDir);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k = 1:nCases
            [caseResults{k},caseFiles{k}] = run_gap_case( ...
                MVec(taskM(k)),gapVec(taskG(k)),taskR(k),o,runMeta,caseDir);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    caseGrid = reshape(caseResults,nM,nGap,nRep);
    fileGrid = reshape(caseFiles,nM,nGap,nRep);
    summary  = build_summary(caseGrid,MVec,gapVec,nRep);

    bundle = struct();
    bundle.meta       = runMeta;
    bundle.experiment = 'minimum_spacing_sensitivity';
    bundle.run = struct( ...
        'id',runId, ...
        'label',char(string(o.RunLabel)), ...
        'conditionTag',conditionTag, ...
        'outputDirectory',outDir);
    bundle.axes = struct('M',MVec,'minGap',gapVec,'replicate',1:nRep, ...
        'SNR_dB',double(o.SNR_dB),'sigma_X_sq',double(o.sigma_X_sq),'P_avg',double(o.P_avg));
    bundle.settings = struct( ...
        'ghN_h',o.ghN_h,'saMaxIter',o.saMaxIter,'saNStarts',o.saNStarts, ...
        'nReplicates',nRep,'baseSeed',o.baseSeed,'historyEvery',o.historyEvery, ...
        'pinZero',o.pinZero,'pairedSeedsAcrossMinGap',true, ...
        'computeBER',o.ComputeBER,'runLabel',char(string(o.RunLabel)));
    bundle.caseFiles = fileGrid;
    bundle.summary   = summary;
    bundle.figures   = {};

    if o.MakePlots
        bundle.figures = make_summary_plots(bundle,figDir,logical(o.SaveFigures));
    end

    bundlePath = fullfile(outDir,'minGapSensitivity.mat');
    save(bundlePath,'bundle','-v7.3');

    fprintf('\nSaved summary bundle: %s\n',bundlePath);
    print_summary(bundle);

    function progress_tick()
        doneCount = doneCount + 1;
        fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function [c,path] = run_gap_case(M,minGap,repIdx,o,runMeta,caseDir)
    % Paired design: seed depends on M/channel/replicate, deliberately NOT d_min.
    replicateBaseSeed = double(o.baseSeed) + (repIdx-1)*10000019;
    seed = fso_result_utils.case_seed(replicateBaseSeed,M,o.SNR_dB,o.sigma_X_sq);

    xMaxBound = fso_result_utils.feasible_xmax_bound(M,o.P_avg,minGap,o.pinZero);
    cfg = build_fso_config(M,o.P_avg,o.SNR_dB,o.sigma_X_sq, ...
        'ghN_h',o.ghN_h, ...
        'xMaxBound',xMaxBound, ...
        'saMaxIter',o.saMaxIter, ...
        'saNStarts',o.saNStarts, ...
        'saUseParallel',false, ...
        'minGap',minGap, ...
        'pinZero',o.pinZero, ...
        'seedInit',seed, ...
        'logEvery',0, ...
        'historyEvery',o.historyEvery);

    xPAM = define_constellation(M,o.P_avg,true,"mean");
    AMI_functions.assert_constellation_feasible(xPAM,cfg,1e-10);

    amiPAMFast = cfg.AMI_Evaluator(xPAM);
    tBase = tic;
    amiPAMVal = cfg.AMI_Validator(xPAM);
    baseValRuntime = toc(tBase);

    label = sprintf('[gap|M%d|d%.4g|rep%d]',M,minGap,repIdx);
    [saOut,starts] = sa_multistart(cfg,xPAM,label);

    xGS       = saOut.bestXValidated(:);
    amiGSVal  = saOut.bestMIValidated;
    amiGSFast = saOut.validatedWinnerFastMI;

    if numel(xGS) > 1
        actualMinGap = min(diff(xGS));
    else
        actualMinGap = NaN;
    end

    startVals = [starts.bestMIValidated];
    restartMean = mean(startVals,'omitnan');
    restartStd  = std(startVals,0,'omitnan');
    restartRange = max(startVals)-min(startVals);

    winStart = starts(saOut.bestStartValidated);
    convIter99 = convergence_iteration(winStart.history,winStart.initialMI,winStart.bestMI);

    skeletonMean = (M-1)*minGap/2;
    constraintLoadFraction = skeletonMean/o.P_avg;

    berPAM = NaN;
    berGS  = NaN;
    if o.ComputeBER
        berPAM = BER_functions.calculate_BER_noCSI_ML(xPAM(:),cfg);
        berGS  = BER_functions.calculate_BER_noCSI_ML(xGS(:),cfg);
    end

    saOutSaved = saOut;
    if isfield(saOutSaved,'bestRun'), saOutSaved = rmfield(saOutSaved,'bestRun'); end

    c = struct();
    c.meta = runMeta;
    c.case = struct('M',M,'P_avg',o.P_avg,'SNR_dB',o.SNR_dB, ...
        'sigma_X_sq',o.sigma_X_sq,'minGap',minGap,'replicate',repIdx,'seed',seed);
    c.config = fso_result_utils.config_snapshot(cfg);
    c.constraint = struct( ...
        'imposedMinGap',minGap, ...
        'skeletonMeanPower',skeletonMean, ...
        'loadFractionOfMeanPower',constraintLoadFraction, ...
        'xMaxFeasibleBound',xMaxBound, ...
        'fastYGridPoints',numel(cfg.y_grid));
    c.baseline = struct('x',xPAM(:),'amiFast',amiPAMFast, ...
        'amiValidated',amiPAMVal,'ber',berPAM);
    c.starts = starts;
    c.optimizer = saOutSaved;
    c.winner = struct( ...
        'startIndex',saOut.bestStartValidated, ...
        'x',xGS, ...
        'amiFast',amiGSFast, ...
        'amiValidated',amiGSVal, ...
        'actualMinGap',actualMinGap, ...
        'xMax',max(xGS), ...
        'acceptanceRate',winStart.acceptanceRate, ...
        'convergenceIter99',convIter99, ...
        'ber',berGS);
    c.selection = struct( ...
        'fastWinnerIndex',saOut.bestStartFast, ...
        'validatedWinnerIndex',saOut.bestStartValidated, ...
        'changed',saOut.selectionChanged, ...
        'validatedAdvantageOverFastWinner',saOut.validatedAdvantageOverFastWinner, ...
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap, ...
        'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    c.repeatability = struct( ...
        'restartValidatedMean',restartMean, ...
        'restartValidatedStd',restartStd, ...
        'restartValidatedRange',restartRange, ...
        'meanAcceptanceRate',saOut.meanAcceptanceRate);
    c.metrics = struct('gainBits',amiGSVal-amiPAMVal, ...
        'fastGainBits',amiGSFast-amiPAMFast);
    c.runtime = struct('sa',saOut.totalSARuntime, ...
        'validation',saOut.validationTotalRuntime+baseValRuntime, ...
        'total',saOut.totalRuntime+baseValRuntime);

    path = fullfile(caseDir,sprintf('case_M%d_gap%s_rep%02d.mat', ...
        M,num_token(minGap,4),repIdx));
    save(path,'c','-v7.3');
end


function iter99 = convergence_iteration(history,initialMI,bestMI)
    improvement = bestMI-initialMI;
    if ~isfinite(improvement) || improvement <= 1e-12
        iter99 = 0;
        return;
    end
    target = initialMI + 0.99*improvement;
    idx = find(history.bestMI >= target,1,'first');
    if isempty(idx)
        iter99 = history.iter(end);
    else
        iter99 = history.iter(idx);
    end
end


function S = build_summary(caseGrid,MVec,gapVec,nRep)
    nM = numel(MVec); nG = numel(gapVec);
    sz = [nM nG nRep];

    S.amiPAMValidated = nan(sz);
    S.amiGSValidated  = nan(sz);
    S.gainBits        = nan(sz);
    S.actualMinGap    = nan(sz);
    S.xMax            = nan(sz);
    S.constraintLoadFraction = nan(sz);
    S.xMaxFeasibleBound = nan(sz);
    S.fastYGridPoints = nan(sz);
    S.restartValidatedMean  = nan(sz);
    S.restartValidatedStd   = nan(sz);
    S.restartValidatedRange = nan(sz);
    S.winnerAcceptanceRate  = nan(sz);
    S.winnerConvergenceIter99 = nan(sz);
    S.selectionChanged = false(sz);
    S.meanAbsFastValidationGap = nan(sz);
    S.maxAbsFastValidationGap  = nan(sz);
    S.totalRuntime = nan(sz);
    S.xWinner = cell(sz);
    S.seed = nan(sz);

    for iM = 1:nM
        for iG = 1:nG
            for iR = 1:nRep
                c = caseGrid{iM,iG,iR};
                S.amiPAMValidated(iM,iG,iR) = c.baseline.amiValidated;
                S.amiGSValidated(iM,iG,iR)  = c.winner.amiValidated;
                S.gainBits(iM,iG,iR)        = c.metrics.gainBits;
                S.actualMinGap(iM,iG,iR)    = c.winner.actualMinGap;
                S.xMax(iM,iG,iR)            = c.winner.xMax;
                S.constraintLoadFraction(iM,iG,iR) = c.constraint.loadFractionOfMeanPower;
                S.xMaxFeasibleBound(iM,iG,iR) = c.constraint.xMaxFeasibleBound;
                S.fastYGridPoints(iM,iG,iR) = c.constraint.fastYGridPoints;
                S.restartValidatedMean(iM,iG,iR)  = c.repeatability.restartValidatedMean;
                S.restartValidatedStd(iM,iG,iR)   = c.repeatability.restartValidatedStd;
                S.restartValidatedRange(iM,iG,iR) = c.repeatability.restartValidatedRange;
                S.winnerAcceptanceRate(iM,iG,iR) = c.winner.acceptanceRate;
                S.winnerConvergenceIter99(iM,iG,iR) = c.winner.convergenceIter99;
                S.selectionChanged(iM,iG,iR) = c.selection.changed;
                S.meanAbsFastValidationGap(iM,iG,iR) = c.selection.meanAbsFastValidationGap;
                S.maxAbsFastValidationGap(iM,iG,iR)  = c.selection.maxAbsFastValidationGap;
                S.totalRuntime(iM,iG,iR) = c.runtime.total;
                S.xWinner{iM,iG,iR} = c.winner.x;
                S.seed(iM,iG,iR) = c.case.seed;
            end
        end
    end

    fields = {'amiPAMValidated','amiGSValidated','gainBits','actualMinGap','xMax', ...
        'constraintLoadFraction','xMaxFeasibleBound','fastYGridPoints', ...
        'restartValidatedMean','restartValidatedStd','restartValidatedRange', ...
        'winnerAcceptanceRate','winnerConvergenceIter99', ...
        'meanAbsFastValidationGap','maxAbsFastValidationGap','totalRuntime'};

    S.mean = struct();
    S.stdAcrossReplicates = struct();
    for k = 1:numel(fields)
        f = fields{k};
        S.mean.(f) = mean(S.(f),3,'omitnan');
        S.stdAcrossReplicates.(f) = std(S.(f),0,3,'omitnan');
    end
    S.selectionChangedFraction = mean(double(S.selectionChanged),3);
end


function files = make_summary_plots(bundle,figDir,saveFigures)
    MVec   = bundle.axes.M;
    gapVec = bundle.axes.minGap;
    S      = bundle.summary;
    nM     = numel(MVec);
    files  = {};

    % Figure 1: primary reviewer-response metrics.
    f1 = figure('Name','Minimum-spacing sensitivity: AMI and shaping gain','Color','w');
    tl = tiledlayout(f1,1,2,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
    for iM = 1:nM
        y = reshape(S.mean.amiGSValidated(iM,:),1,[]);
        e = reshape(S.stdAcrossReplicates.amiGSValidated(iM,:),1,[]);
        errorbar(ax1,gapVec,y,e,'-o','LineWidth',1.4,'DisplayName',sprintf('M=%d',MVec(iM)));
    end
    xlabel(ax1,'Imposed minimum spacing d_{min}');
    ylabel(ax1,'Validated GS AMI [bits/symbol]');
    title(ax1,sprintf('SNR=%.1f dB, \\sigma_X^2=%.3f',bundle.axes.SNR_dB,bundle.axes.sigma_X_sq));
    legend(ax1,'Location','best');

    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    for iM = 1:nM
        y = reshape(S.mean.gainBits(iM,:),1,[]);
        e = reshape(S.stdAcrossReplicates.gainBits(iM,:),1,[]);
        errorbar(ax2,gapVec,y,e,'-o','LineWidth',1.4,'DisplayName',sprintf('M=%d',MVec(iM)));
    end
    xlabel(ax2,'Imposed minimum spacing d_{min}');
    ylabel(ax2,'AMI shaping gain [bits/symbol]');
    title(ax2,'Validated GS - Uniform PAM');
    legend(ax2,'Location','best');

    if saveFigures
        files = [files save_figure_pair(f1,figDir,'minGap_AMI_and_gain')]; %#ok<AGROW>
    end

    % Figure 2: whether the constraint is active and how geometry changes.
    f2 = figure('Name','Minimum-spacing sensitivity: geometry diagnostics','Color','w');
    tl2 = tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');

    ax3 = nexttile(tl2,1); hold(ax3,'on'); grid(ax3,'on');
    plot(ax3,gapVec,gapVec,'--','LineWidth',1.0,'DisplayName','actual = imposed');
    for iM = 1:nM
        y = reshape(S.mean.actualMinGap(iM,:),1,[]);
        plot(ax3,gapVec,y,'-o','LineWidth',1.4,'DisplayName',sprintf('M=%d',MVec(iM)));
    end
    xlabel(ax3,'Imposed d_{min}');
    ylabel(ax3,'Actual minimum adjacent spacing');
    title(ax3,'Is the spacing constraint active?');
    legend(ax3,'Location','best');

    ax4 = nexttile(tl2,2); hold(ax4,'on'); grid(ax4,'on');
    for iM = 1:nM
        y = reshape(S.mean.xMax(iM,:),1,[]);
        plot(ax4,gapVec,y,'-o','LineWidth',1.4,'DisplayName',sprintf('M=%d',MVec(iM)));
    end
    xlabel(ax4,'Imposed d_{min}');
    ylabel(ax4,'Largest optimized level x_{max}');
    title(ax4,'Peak level induced by optimized geometry');
    legend(ax4,'Location','best');

    if saveFigures
        files = [files save_figure_pair(f2,figDir,'minGap_geometry_diagnostics')]; %#ok<AGROW>
    end

    % One constellation-geometry figure per M. For multiple replicates use
    % the replicate whose validated AMI is closest to the replicate mean.
    for iM = 1:nM
        f = figure('Name',sprintf('Minimum-spacing geometry M=%d',MVec(iM)),'Color','w');
        ax = axes(f); hold(ax,'on'); grid(ax,'on');
        for iG = 1:numel(gapVec)
            vals = reshape(S.amiGSValidated(iM,iG,:),1,[]);
            target = mean(vals,'omitnan');
            [~,rep] = min(abs(vals-target));
            x = S.xWinner{iM,iG,rep};
            plot(ax,1:numel(x),x,'-o','LineWidth',1.2, ...
                'DisplayName',sprintf('d_{min}=%.4g',gapVec(iG)));
        end
        xlabel(ax,'Ordered constellation level index');
        ylabel(ax,'Optical intensity level');
        title(ax,sprintf('Optimized constellations: M=%d, SNR=%.1f dB, \\sigma_X^2=%.3f', ...
            MVec(iM),bundle.axes.SNR_dB,bundle.axes.sigma_X_sq));
        legend(ax,'Location','best');
        if saveFigures
            files = [files save_figure_pair(f,figDir,sprintf('minGap_constellations_M%d',MVec(iM)))]; %#ok<AGROW>
        end
    end
end


function files = save_figure_pair(fig,figDir,baseName)
    pngPath = fullfile(figDir,[baseName '.png']);
    figPath = fullfile(figDir,[baseName '.fig']);
    exportgraphics(fig,pngPath,'Resolution',200);
    savefig(fig,figPath);
    files = {pngPath,figPath};
end


function print_summary(bundle)
    MVec = bundle.axes.M;
    gapVec = bundle.axes.minGap;
    S = bundle.summary;

    fprintf('\nValidated minimum-spacing sensitivity summary\n');
    fprintf('     M |    d_min | GS AMI mean | gain mean | actual gap | x_max | constraint load | restart std\n');
    fprintf('-----------------------------------------------------------------------------------------------\n');
    for iM = 1:numel(MVec)
        for iG = 1:numel(gapVec)
            fprintf('%6d | %8.4f | %11.6f | %+9.6f | %10.5f | %5.2f | %14.3f | %.3e\n', ...
                MVec(iM),gapVec(iG), ...
                S.mean.amiGSValidated(iM,iG), ...
                S.mean.gainBits(iM,iG), ...
                S.mean.actualMinGap(iM,iG), ...
                S.mean.xMax(iM,iG), ...
                S.mean.constraintLoadFraction(iM,iG), ...
                S.mean.restartValidatedStd(iM,iG));
        end
    end
    fprintf('\nselectionChanged fraction across all cases: %.3f\n',mean(S.selectionChanged(:)));
    fprintf('Max fast-vs-validation gap: %.3e bits/symbol\n', ...
        max(S.maxAbsFastValidationGap(:),[],'omitnan'));
end


function token = sanitize_run_label(v)
    token = strtrim(char(string(v)));
    if isempty(token)
        token = 'run';
        return;
    end
    token = regexprep(token,'[^A-Za-z0-9._-]+','_');
    token = regexprep(token,'_+','_');
    token = regexprep(token,'^[._-]+|[._-]+$','');
    if isempty(token), token = 'run'; end
end


function [outDir,runId] = make_unique_run_directory(conditionDir,requestedRunId)
    if ~exist(conditionDir,'dir'), mkdir(conditionDir); end

    runId = requestedRunId;
    outDir = fullfile(conditionDir,runId);
    suffix = 1;
    while exist(outDir,'dir')
        suffix = suffix + 1;
        runId = sprintf('%s_%02d',requestedRunId,suffix);
        outDir = fullfile(conditionDir,runId);
    end

    [ok,msg] = mkdir(outDir);
    if ~ok
        error('sim_min_gap_sensitivity:RunDirectoryCreateFailed', ...
            'Could not create isolated run directory %s: %s',outDir,msg);
    end
end


function token = num_token(v,nDigits)
    if nargin < 2, nDigits = 4; end
    fmt = sprintf('%%+.%df',nDigits);
    token = sprintf(fmt,v);
    token = strrep(token,'+','p');
    token = strrep(token,'-','m');
    token = strrep(token,'.','p');
end