function dbundle = diagnose_y_grid_runtime(source,varargin)
%DIAGNOSE_Y_GRID_RUNTIME  Commit I.5 candidate-adaptive y-grid diagnostic.
%
% This experiment is intentionally diagnostic-only.  It does NOT modify the
% production evaluator, the SA algorithm, Commit G/I.4, or any constellation.
% It asks whether the fixed fast-AMI y-grid (currently built once from the
% global feasible xMax bound) can be replaced later by a candidate-adaptive
% y-grid built from the actual maximum level of the constellation being
% scored.
%
% For every selected validated winner from sim_revision_AMI_vs_SNR(), the
% function compares:
%
%   REFERENCE (current production path)
%       yGrid = build_noCSI_y_grid(cfg, feasibleXMaxBound)
%
%   ADAPTIVE candidate path
%       yGrid = build_noCSI_y_grid(cfg, scale * max(abs(x)))
%
% for several support scale factors.  The adaptive timing INCLUDES y-grid
% construction on every call, so the reported speedup is conservative and
% directly relevant to a future SA objective implementation.
%
% IMPORTANT: scale*max(x) is used only to choose numerical integration
% support.  It does not project, clip, reject, or constrain x, and therefore
% does NOT introduce a peak-power constraint.  Every candidate is scored on
% a grid tailored to that exact candidate.
%
% The log-h fading quadrature, dt, t-span, outer-grid construction formula,
% and independent validator are otherwise unchanged.
%
% SOURCE may be either:
%   * an in-memory bundle from sim_revision_AMI_vs_SNR(), or
%   * a path to revisedAMIvsSNR.mat.
%
% Name-value options:
%   'MVec'              M values ([] -> all source values)
%   'SNRVec'            SNR values ([] -> all source values)
%   'TurbulenceVec'     sigma_X^2 values ([] -> all source values)
%   'ReplicateVec'      replicate indices ([] -> all source replicates)
%   'ScaleVec'          actual-xMax support multipliers
%                       (default [1 1.1 1.25 1.5])
%   'AbsErrorTolerance' all-case target versus validator (default 1e-3)
%   'TimingWarmups'     untimed warmups per method/candidate (default 1)
%   'TimingRepeats'     timed evaluations per method/candidate (default 3)
%   'loghDt'            [] -> inherit source method (default [])
%   'loghSpanSigma'     [] -> inherit source method (default [])
%   'loghYBlockSize'    [] -> inherit source method (default [])
%   'UseParallelCases'  parallelize candidate cases (default true)
%   'MakePlots'         create summary figure (default true)
%   'SaveFigures'       save PNG/FIG pair (default true)
%   'RunLabel'          optional run label
%   'ResultsRoot'       output root
%
% Recommended first diagnostic on the Commit-I.3 pilot:
%   d5 = diagnose_y_grid_runtime(i3Pilot, ...
%       'MVec',[8 32], 'SNRVec',[5 20 30], ...
%       'TurbulenceVec',0.1, 'ReplicateVec',[1 2], ...
%       'ScaleVec',[1 1.1 1.25 1.5], ...
%       'TimingWarmups',1, 'TimingRepeats',3, ...
%       'UseParallelCases',true, 'RunLabel','I5_grid_pilot');
%
% A future production change should only be considered if an adaptive scale
% preserves the desired fast-vs-validator accuracy on every selected case
% while materially reducing total evaluation runtime.

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p,'source',@(v) isstruct(v) || ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'MVec',[],@(v) isempty(v) || numeric_vector(v));
    addParameter(p,'SNRVec',[],@(v) isempty(v) || numeric_vector(v));
    addParameter(p,'TurbulenceVec',[],@(v) isempty(v) || (numeric_vector(v) && all(v>=0)));
    addParameter(p,'ReplicateVec',[],@(v) isempty(v) || (numeric_vector(v) && all(v>=1) && all(mod(v,1)==0)));
    addParameter(p,'ScaleVec',[1 1.1 1.25 1.5],@(v) numeric_vector(v) && all(v>=1));
    addParameter(p,'AbsErrorTolerance',1e-3,@positive_scalar);
    addParameter(p,'TimingWarmups',1,@nonnegative_integer);
    addParameter(p,'TimingRepeats',3,@positive_integer);
    addParameter(p,'loghDt',[],@(v) isempty(v) || positive_scalar(v));
    addParameter(p,'loghSpanSigma',[],@(v) isempty(v) || positive_scalar(v));
    addParameter(p,'loghYBlockSize',[],@(v) isempty(v) || positive_integer(v));
    addParameter(p,'UseParallelCases',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'MakePlots',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'SaveFigures',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v) || isstring(v));
    parse(p,source,varargin{:});
    o = p.Results;

    [src,sourceDescription] = load_source(source);
    validate_source(src);

    srcM = double(src.axes.M(:).');
    srcSNR = double(src.axes.SNR_dB(:).');
    srcSig = double(src.axes.sigma_X_sq(:).');
    srcRep = double(src.axes.replicate(:).');

    if isempty(o.MVec), MVec=srcM; else, MVec=unique(double(o.MVec(:).'),'stable'); end
    if isempty(o.SNRVec), snrVec=srcSNR; else, snrVec=unique(double(o.SNRVec(:).'),'stable'); end
    if isempty(o.TurbulenceVec), sigVec=srcSig; else, sigVec=unique(double(o.TurbulenceVec(:).'),'stable'); end
    if isempty(o.ReplicateVec), repVec=srcRep; else, repVec=unique(double(o.ReplicateVec(:).'),'stable'); end
    scaleVec = unique(double(o.ScaleVec(:).'),'stable');

    idxM = map_axis(MVec,srcM,'M');
    idxSNR = map_axis(snrVec,srcSNR,'SNR');
    idxSig = map_axis(sigVec,srcSig,'sigma_X_sq');
    idxRep = map_axis(repVec,srcRep,'replicate');

    [dt,spanSigma,yBlockSize] = resolve_logh_settings(src,o);

    nM=numel(MVec); nS=numel(snrVec); nG=numel(sigVec); nR=numel(repVec);
    [iMGrid,iGGrid,iSGrid,iRGrid] = ndgrid(1:nM,1:nG,1:nS,1:nR);
    taskM=iMGrid(:); taskG=iGGrid(:); taskS=iSGrid(:); taskR=iRGrid(:);
    nCases=numel(taskM);

    % Immutable result directory.
    resultsRoot=char(o.ResultsRoot);
    parentDir=fullfile(resultsRoot,'y_grid_runtime_diagnostic');
    labelToken=sanitize_label(o.RunLabel);
    scaleToken=['S' join_numbers(scaleVec,2)];
    methodToken=sprintf('LOGHdt%s_T%s',num_token(dt,4),num_token(spanSigma,2));
    utcToken=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requestedRunId=sprintf('%s_%s_%s_%s',labelToken,methodToken,scaleToken,utcToken);
    [outDir,runId]=make_unique_run_directory(parentDir,requestedRunId);
    caseDir=fullfile(outDir,'cases'); mkdir(caseDir);
    figDir=fullfile(outDir,'figures');
    if o.MakePlots && o.SaveFigures, mkdir(figDir); end

    runMeta=fso_result_utils.run_metadata(mfilename);
    runMeta.runId=runId;
    runMeta.runLabel=char(string(o.RunLabel));
    runMeta.outputDirectory=outDir;

    fprintf('\n============================================================\n');
    fprintf('Commit I.5 - candidate-adaptive y-grid runtime diagnostic\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('Source: %s\n',sourceDescription);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | reps=%s\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),mat2str(repVec));
    fprintf('log-h: dt=%.5g | span=+/-%.3g sigma_t | yBlock=%d\n',dt,spanSigma,yBlockSize);
    fprintf('Adaptive xMax support scales=%s\n',mat2str(scaleVec));
    fprintf('Accuracy target: max |fast-validator| <= %.3e bit/symbol\n',o.AbsErrorTolerance);
    fprintf('Timing: warmups=%d | repeats=%d | adaptive timing includes grid construction\n', ...
        o.TimingWarmups,o.TimingRepeats);
    fprintf('Cases=%d | adaptive evaluations=%d (+ fixed references)\n', ...
        nCases,nCases*numel(scaleVec));
    fprintf('No constellation clipping/projection/peak constraint is introduced.\n');
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    caseResults=cell(nCases,1);
    caseFiles=cell(nCases,1);

    doPar=logical(o.UseParallelCases);
    if doPar
        try
            doPar=license('test','Distrib_Computing_Toolbox')~=0;
        catch
            doPar=false;
        end
    end

    if doPar
        pool=gcp('nocreate');
        if isempty(pool), pool=parpool('Processes'); end
        fprintf('Outer case parallelism: %d process workers\n',pool.NumWorkers);
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~) progress_tick());
        parfor k=1:nCases
            [caseResults{k},caseFiles{k}] = run_case( ...
                taskM(k),taskG(k),taskS(k),taskR(k), ...
                MVec,sigVec,snrVec,repVec,idxM,idxSig,idxSNR,idxRep, ...
                scaleVec,src,dt,spanSigma,yBlockSize, ...
                double(o.TimingWarmups),double(o.TimingRepeats),runMeta,caseDir);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k=1:nCases
            [caseResults{k},caseFiles{k}] = run_case( ...
                taskM(k),taskG(k),taskS(k),taskR(k), ...
                MVec,sigVec,snrVec,repVec,idxM,idxSig,idxSNR,idxRep, ...
                scaleVec,src,dt,spanSigma,yBlockSize, ...
                double(o.TimingWarmups),double(o.TimingRepeats),runMeta,caseDir);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    caseGrid=reshape(caseResults,nM,nG,nS,nR);
    fileGrid=reshape(caseFiles,nM,nG,nS,nR);
    summary=build_summary(caseGrid,scaleVec,double(o.AbsErrorTolerance));

    dbundle=struct();
    dbundle.meta=runMeta;
    dbundle.experiment='candidate_adaptive_y_grid_runtime_diagnostic';
    dbundle.run=struct('id',runId,'label',char(string(o.RunLabel)), ...
        'outputDirectory',outDir);
    dbundle.source=struct('description',sourceDescription,'runId',source_run_id(src));
    dbundle.axes=struct('M',MVec,'SNR_dB',snrVec,'sigma_X_sq',sigVec, ...
        'replicate',repVec,'supportScale',scaleVec);
    dbundle.method=struct('fastFadingMethod','logh','loghDt',dt, ...
        'loghSpanSigma',spanSigma,'loghYBlockSize',yBlockSize, ...
        'absErrorTolerance',double(o.AbsErrorTolerance), ...
        'timingWarmups',double(o.TimingWarmups),'timingRepeats',double(o.TimingRepeats));
    dbundle.caseFiles=fileGrid;
    dbundle.summary=summary;
    dbundle.figureFiles={};

    if o.MakePlots
        dbundle.figureFiles=make_plots(dbundle,figDir,logical(o.SaveFigures));
    end

    bundlePath=fullfile(outDir,'yGridRuntimeDiagnostic.mat');
    save(bundlePath,'dbundle','-v7.3');
    fprintf('\nSaved y-grid diagnostic bundle: %s\n',bundlePath);
    print_summary(dbundle,caseGrid);

    function progress_tick()
        doneCount=doneCount+1;
        fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function [c,path] = run_case(iM,iG,iS,iR,MVec,sigVec,snrVec,repVec, ...
        idxM,idxSig,idxSNR,idxRep,scaleVec,src,dt,spanSigma,yBlockSize, ...
        nWarm,nRepeat,runMeta,caseDir)

    M=MVec(iM); sig=sigVec(iG); snr=snrVec(iS); rep=repVec(iR);
    sM=idxM(iM); sG=idxSig(iG); sS=idxSNR(iS); sR=idxRep(iR);

    x=real(src.summary.xWinner{sM,sG,sS,sR}(:));
    storedVal=src.summary.raw.amiGSValidated(sM,sG,sS,sR);
    Pavg=double(src.axes.P_avg);
    minGap=double(src.method.minGap);
    pinZero=logical(src.method.pinZero);

    if isfield(src,'xMaxBounds') && numel(src.xMaxBounds)>=sM
        xMaxBound=double(src.xMaxBounds(sM));
    else
        xMaxBound=fso_result_utils.feasible_xmax_bound(M,Pavg,minGap,pinZero);
    end

    cfg=build_fso_config(M,Pavg,snr,sig, ...
        'FastFadingMethod','logh','loghDt',dt,'loghSpanSigma',spanSigma, ...
        'loghYBlockSize',yBlockSize,'xMaxBound',xMaxBound, ...
        'saMaxIter',1,'saNStarts',1,'saUseParallel',false, ...
        'minGap',minGap,'pinZero',pinZero,'seedInit',1,'logEvery',0,'historyEvery',1);

    AMI_functions.assert_constellation_feasible(x,cfg,1e-10);
    validatorAMI=cfg.AMI_Validator(x);
    validatorMismatch=validatorAMI-storedVal;

    fixedGrid=cfg.y_grid;
    [fixedAMI,fixedRuntime] = time_fixed(x,cfg,dt,spanSigma,fixedGrid,yBlockSize,nWarm,nRepeat);
    fixedError=fixedAMI-validatorAMI;

    nScale=numel(scaleVec);
    adaptiveAMI=nan(1,nScale); adaptiveError=nan(1,nScale);
    adaptiveRuntime=nan(1,nScale); adaptiveBuildRuntime=nan(1,nScale);
    adaptiveEvalRuntime=nan(1,nScale); adaptiveNy=nan(1,nScale);
    adaptiveYMin=nan(1,nScale); adaptiveYMax=nan(1,nScale);
    supportXMax=nan(1,nScale);

    actualXMax=max(abs(x));
    for j=1:nScale
        supportXMax(j)=scaleVec(j)*actualXMax;
        [adaptiveAMI(j),adaptiveRuntime(j),adaptiveBuildRuntime(j), ...
            adaptiveEvalRuntime(j),gridInfo] = time_adaptive( ...
            x,cfg,dt,spanSigma,yBlockSize,supportXMax(j),nWarm,nRepeat);
        adaptiveError(j)=adaptiveAMI(j)-validatorAMI;
        adaptiveNy(j)=gridInfo.N;
        adaptiveYMin(j)=gridInfo.yMin;
        adaptiveYMax(j)=gridInfo.yMax;
    end

    c=struct();
    c.meta=runMeta;
    c.case=struct('M',M,'SNR_dB',snr,'sigma_X_sq',sig,'replicate',rep);
    c.source=struct('x',x,'storedValidatedAMI',storedVal, ...
        'recomputedValidatorAMI',validatorAMI,'validatorMismatch',validatorMismatch, ...
        'actualXMax',actualXMax,'feasibleXMaxBound',xMaxBound, ...
        'xMaxToBoundRatio',actualXMax/xMaxBound);
    c.reference=struct('ami',fixedAMI,'signedError',fixedError,'absError',abs(fixedError), ...
        'runtimeSeconds',fixedRuntime,'gridN',numel(fixedGrid), ...
        'gridMin',fixedGrid(1),'gridMax',fixedGrid(end));
    c.adaptive=struct('scale',scaleVec,'supportXMax',supportXMax, ...
        'ami',adaptiveAMI,'signedError',adaptiveError,'absError',abs(adaptiveError), ...
        'runtimeSeconds',adaptiveRuntime,'gridBuildSeconds',adaptiveBuildRuntime, ...
        'evaluatorSeconds',adaptiveEvalRuntime,'gridN',adaptiveNy, ...
        'gridMin',adaptiveYMin,'gridMax',adaptiveYMax, ...
        'speedupVsFixed',fixedRuntime./adaptiveRuntime, ...
        'gridNRatioVsFixed',adaptiveNy/numel(fixedGrid));

    fname=sprintf('ygrid_M%d_SNR%s_sig%s_rep%02d.mat', ...
        M,num_token(snr,2),num_token(sig,4),rep);
    path=fullfile(caseDir,fname);
    save(path,'c','-v7.3');
end


function [ami,medianRuntime] = time_fixed(x,cfg,dt,spanSigma,yGrid,yBlockSize,nWarm,nRepeat)
    for k=1:nWarm
        AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,dt,spanSigma,yGrid,yBlockSize);
    end
    times=nan(1,nRepeat); vals=nan(1,nRepeat);
    for k=1:nRepeat
        t=tic;
        vals(k)=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,dt,spanSigma,yGrid,yBlockSize);
        times(k)=toc(t);
    end
    ami=median(vals,'omitnan');
    medianRuntime=median(times,'omitnan');
end


function [ami,medianTotal,medianBuild,medianEval,gridInfo] = time_adaptive( ...
        x,cfg,dt,spanSigma,yBlockSize,supportXMax,nWarm,nRepeat)

    for k=1:nWarm
        yg=AMI_functions.build_noCSI_y_grid(cfg,supportXMax);
        AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,dt,spanSigma,yg,yBlockSize);
    end

    totalTimes=nan(1,nRepeat); buildTimes=nan(1,nRepeat); evalTimes=nan(1,nRepeat);
    vals=nan(1,nRepeat); lastGrid=[];
    for k=1:nRepeat
        tAll=tic;
        tBuild=tic;
        yg=AMI_functions.build_noCSI_y_grid(cfg,supportXMax);
        buildTimes(k)=toc(tBuild);
        tEval=tic;
        vals(k)=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,dt,spanSigma,yg,yBlockSize);
        evalTimes(k)=toc(tEval);
        totalTimes(k)=toc(tAll);
        lastGrid=yg;
    end

    ami=median(vals,'omitnan');
    medianTotal=median(totalTimes,'omitnan');
    medianBuild=median(buildTimes,'omitnan');
    medianEval=median(evalTimes,'omitnan');
    gridInfo=struct('N',numel(lastGrid),'yMin',lastGrid(1),'yMax',lastGrid(end));
end


function S = build_summary(caseGrid,scaleVec,tol)
    nCases=numel(caseGrid); nScale=numel(scaleVec);
    fixedAbs=nan(1,nCases); fixedSigned=nan(1,nCases); fixedRun=nan(1,nCases);
    fixedN=nan(1,nCases); ratios=nan(1,nCases); valMismatch=nan(1,nCases);
    adaptAbs=nan(nScale,nCases); adaptSigned=nan(nScale,nCases);
    adaptRun=nan(nScale,nCases); adaptBuild=nan(nScale,nCases);
    adaptEval=nan(nScale,nCases); adaptN=nan(nScale,nCases);
    speedup=nan(nScale,nCases); nRatio=nan(nScale,nCases);

    for k=1:nCases
        c=caseGrid{k};
        fixedAbs(k)=c.reference.absError;
        fixedSigned(k)=c.reference.signedError;
        fixedRun(k)=c.reference.runtimeSeconds;
        fixedN(k)=c.reference.gridN;
        ratios(k)=c.source.xMaxToBoundRatio;
        valMismatch(k)=c.source.validatorMismatch;
        adaptAbs(:,k)=c.adaptive.absError(:);
        adaptSigned(:,k)=c.adaptive.signedError(:);
        adaptRun(:,k)=c.adaptive.runtimeSeconds(:);
        adaptBuild(:,k)=c.adaptive.gridBuildSeconds(:);
        adaptEval(:,k)=c.adaptive.evaluatorSeconds(:);
        adaptN(:,k)=c.adaptive.gridN(:);
        speedup(:,k)=c.adaptive.speedupVsFixed(:);
        nRatio(:,k)=c.adaptive.gridNRatioVsFixed(:);
    end

    S=struct();
    S.reference=struct( ...
        'meanAbsError',mean(fixedAbs,'omitnan'), ...
        'maxAbsError',max(fixedAbs,[],'omitnan'), ...
        'meanSignedError',mean(fixedSigned,'omitnan'), ...
        'medianRuntimeSeconds',median(fixedRun,'omitnan'), ...
        'meanRuntimeSeconds',mean(fixedRun,'omitnan'), ...
        'meanGridN',mean(fixedN,'omitnan'), ...
        'maxGridN',max(fixedN,[],'omitnan'));

    A=struct();
    A.scale=scaleVec;
    A.meanAbsError=mean(adaptAbs,2,'omitnan').';
    A.maxAbsError=max(adaptAbs,[],2,'omitnan').';
    A.meanSignedError=mean(adaptSigned,2,'omitnan').';
    A.medianRuntimeSeconds=median(adaptRun,2,'omitnan').';
    A.meanRuntimeSeconds=mean(adaptRun,2,'omitnan').';
    A.medianGridBuildSeconds=median(adaptBuild,2,'omitnan').';
    A.medianEvaluatorSeconds=median(adaptEval,2,'omitnan').';
    A.meanGridN=mean(adaptN,2,'omitnan').';
    A.maxGridN=max(adaptN,[],2,'omitnan').';
    A.medianSpeedupVsFixed=median(speedup,2,'omitnan').';
    A.meanSpeedupVsFixed=mean(speedup,2,'omitnan').';
    A.minSpeedupVsFixed=min(speedup,[],2,'omitnan').';
    A.meanGridNRatioVsFixed=mean(nRatio,2,'omitnan').';
    A.maxGridNRatioVsFixed=max(nRatio,[],2,'omitnan').';
    A.allCasesPassTolerance=A.maxAbsError<=tol;
    S.adaptive=A;

    pass=find(A.allCasesPassTolerance);
    if isempty(pass)
        S.smallestScaleAllCasesPass=NaN;
        S.fastestPassingScale=NaN;
    else
        [~,ii]=min(scaleVec(pass));
        S.smallestScaleAllCasesPass=scaleVec(pass(ii));
        [~,jj]=max(A.medianSpeedupVsFixed(pass));
        S.fastestPassingScale=scaleVec(pass(jj));
    end

    S.absErrorTolerance=tol;
    S.maxAbsValidatorRecomputeMismatch=max(abs(valMismatch),[],'omitnan');
    S.meanXMaxToBoundRatio=mean(ratios,'omitnan');
    S.minXMaxToBoundRatio=min(ratios,[],'omitnan');
    S.maxXMaxToBoundRatio=max(ratios,[],'omitnan');
end


function files = make_plots(bundle,figDir,saveFigures)
    s=bundle.axes.supportScale;
    A=bundle.summary.adaptive;
    tol=bundle.summary.absErrorTolerance;

    f=figure('Name','Commit I.5 y-grid diagnostic','Color','w');
    tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
    title(tl,'Candidate-adaptive y-grid diagnostic');

    ax1=nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
    semilogy(ax1,s,A.maxAbsError,'-o','LineWidth',1.5);
    yline(ax1,tol,'--','Accuracy target');
    xlabel(ax1,'Support scale x actual x_{max}'); ylabel(ax1,'Max |fast-validator| [bit/symbol]');
    title(ax1,'Accuracy');

    ax2=nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
    plot(ax2,s,A.medianSpeedupVsFixed,'-o','LineWidth',1.5);
    yline(ax2,1,'--','No speedup');
    xlabel(ax2,'Support scale'); ylabel(ax2,'Median speedup vs fixed-bound grid');
    title(ax2,'Runtime');

    ax3=nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
    plot(ax3,s,A.meanGridNRatioVsFixed,'-o','LineWidth',1.5);
    yline(ax3,1,'--','Current grid');
    xlabel(ax3,'Support scale'); ylabel(ax3,'Mean N_y / current N_y');
    title(ax3,'Grid size');

    files={};
    if saveFigures
        if ~exist(figDir,'dir'), mkdir(figDir); end
        pngPath=fullfile(figDir,'candidate_adaptive_y_grid.png');
        figPath=fullfile(figDir,'candidate_adaptive_y_grid.fig');
        exportgraphics(f,pngPath,'Resolution',220);
        savefig(f,figPath);
        files={pngPath,figPath};
    end
end


function print_summary(bundle,caseGrid)
    S=bundle.summary; A=S.adaptive;
    fprintf('\nCANDIDATE-ADAPTIVE Y-GRID SUMMARY\n');
    fprintf('Reference fixed-bound grid:\n');
    fprintf('  mean |error| = %.3e | max |error| = %.3e\n', ...
        S.reference.meanAbsError,S.reference.maxAbsError);
    fprintf('  median runtime = %.4fs | mean N_y = %.1f | max N_y = %.0f\n', ...
        S.reference.medianRuntimeSeconds,S.reference.meanGridN,S.reference.maxGridN);
    fprintf('\nAdaptive grids (timing includes grid construction):\n');
    fprintf('  scale      max|err|      median[s]   med speedup   mean N_y ratio   pass\n');
    for j=1:numel(A.scale)
        fprintf('  %-8.3g   %-12.3e  %-10.4f  %-11.3fx  %-14.3f  %d\n', ...
            A.scale(j),A.maxAbsError(j),A.medianRuntimeSeconds(j), ...
            A.medianSpeedupVsFixed(j),A.meanGridNRatioVsFixed(j),A.allCasesPassTolerance(j));
    end
    fprintf('\nValidator recompute max mismatch: %.3e bit/symbol\n', ...
        S.maxAbsValidatorRecomputeMismatch);
    fprintf('Actual xMax / feasible-bound ratio: mean=%.3f, range=[%.3f %.3f]\n', ...
        S.meanXMaxToBoundRatio,S.minXMaxToBoundRatio,S.maxXMaxToBoundRatio);
    if isnan(S.smallestScaleAllCasesPass)
        fprintf('No tested adaptive scale met the all-case tolerance %.3e.\n',S.absErrorTolerance);
    else
        fprintf('Smallest scale meeting tolerance in every case: %.3g\n',S.smallestScaleAllCasesPass);
        fprintf('Fastest passing scale (median speedup criterion): %.3g\n',S.fastestPassingScale);
    end

    fprintf('\nPer-case reference -> adaptive(scale=1) snapshot:\n');
    for k=1:numel(caseGrid)
        c=caseGrid{k};
        j=find(abs(c.adaptive.scale-1)<1e-12,1,'first');
        if isempty(j), j=1; end
        fprintf('  M=%d SNR=%.1f sig=%.3f rep=%d | xMax/bound=%.3f | N_y %d->%d | err %.2e->%.2e | speedup %.2fx\n', ...
            c.case.M,c.case.SNR_dB,c.case.sigma_X_sq,c.case.replicate, ...
            c.source.xMaxToBoundRatio,c.reference.gridN,c.adaptive.gridN(j), ...
            c.reference.absError,c.adaptive.absError(j),c.adaptive.speedupVsFixed(j));
    end
end


function [dt,spanSigma,yBlockSize] = resolve_logh_settings(src,o)
    if ~isfield(src,'method')
        error('diagnose_y_grid_runtime:MissingMethod','Source bundle has no method metadata.');
    end
    if isfield(src.method,'fastFadingMethod') && ~strcmpi(src.method.fastFadingMethod,'logh')
        error('diagnose_y_grid_runtime:SourceNotLogh', ...
            'Source fastFadingMethod is %s; I.5 expects a log-h production bundle.',src.method.fastFadingMethod);
    end

    if isempty(o.loghDt)
        require_field(src.method,'loghDt'); dt=double(src.method.loghDt);
    else
        dt=double(o.loghDt);
    end
    if isempty(o.loghSpanSigma)
        require_field(src.method,'loghSpanSigma'); spanSigma=double(src.method.loghSpanSigma);
    else
        spanSigma=double(o.loghSpanSigma);
    end
    if isempty(o.loghYBlockSize)
        require_field(src.method,'loghYBlockSize'); yBlockSize=double(src.method.loghYBlockSize);
    else
        yBlockSize=double(o.loghYBlockSize);
    end
end


function [src,description] = load_source(source)
    if isstruct(source)
        src=source; description='<workspace bundle>'; return;
    end
    path=char(source);
    if ~isfile(path)
        error('diagnose_y_grid_runtime:SourceNotFound','Source MAT file not found: %s',path);
    end
    tmp=load(path);
    if isfield(tmp,'bundle')
        src=tmp.bundle;
    elseif isfield(tmp,'i3Pilot')
        src=tmp.i3Pilot;
    else
        names=fieldnames(tmp);
        hit=[];
        for k=1:numel(names)
            if isstruct(tmp.(names{k})) && isfield(tmp.(names{k}),'experiment')
                hit=k; break;
            end
        end
        if isempty(hit)
            error('diagnose_y_grid_runtime:BadSourceFile','No revised AMI-vs-SNR bundle found in %s.',path);
        end
        src=tmp.(names{hit});
    end
    description=path;
end


function validate_source(src)
    required={'axes','method','summary'};
    for k=1:numel(required)
        if ~isfield(src,required{k})
            error('diagnose_y_grid_runtime:BadSource','Source missing field %s.',required{k});
        end
    end
    if ~isfield(src.summary,'xWinner') || ~isfield(src.summary,'raw') || ...
            ~isfield(src.summary.raw,'amiGSValidated')
        error('diagnose_y_grid_runtime:BadSource','Source summary lacks winners/raw validated AMI.');
    end
    axisFields={'M','P_avg','SNR_dB','sigma_X_sq','replicate'};
    for k=1:numel(axisFields)
        if ~isfield(src.axes,axisFields{k})
            error('diagnose_y_grid_runtime:BadSource','Source axes missing %s.',axisFields{k});
        end
    end
    if ~isfield(src.method,'minGap') || ~isfield(src.method,'pinZero')
        error('diagnose_y_grid_runtime:BadSource','Source method lacks minGap/pinZero.');
    end
end


function idx = map_axis(requested,available,label)
    idx=nan(size(requested));
    for k=1:numel(requested)
        tol=1e-10*max(1,abs(requested(k)));
        hit=find(abs(available-requested(k))<=tol,1,'first');
        if isempty(hit)
            error('diagnose_y_grid_runtime:AxisValueNotFound', ...
                'Requested %s value %.12g is not present in source.',label,requested(k));
        end
        idx(k)=hit;
    end
end


function require_field(s,name)
    if ~isfield(s,name)
        error('diagnose_y_grid_runtime:MissingSetting','Source method is missing %s.',name);
    end
end


function id = source_run_id(src)
    id='unknown';
    if isfield(src,'run') && isfield(src.run,'id'), id=src.run.id; end
end


function tf = numeric_vector(v)
    tf=isnumeric(v) && isvector(v) && ~isempty(v) && all(isfinite(v));
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


function token = sanitize_label(label)
    token=strtrim(char(string(label)));
    if isempty(token), token='run'; return; end
    token=regexprep(token,'[^A-Za-z0-9._-]+','_');
    token=regexprep(token,'_+','_');
    token=regexprep(token,'^[._-]+|[._-]+$','');
    if isempty(token), token='run'; end
end


function token = join_numbers(v,nDec)
    parts=cell(1,numel(v));
    for k=1:numel(v), parts{k}=num_token(v(k),nDec); end
    token=strjoin(parts,'-');
end


function token = num_token(v,nDec)
    token=sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token=strrep(token,'-','m');
    token=strrep(token,'+','p');
    token=strrep(token,'.','p');
end


function [outDir,runId] = make_unique_run_directory(parentDir,requestedRunId)
    if ~exist(parentDir,'dir'), mkdir(parentDir); end
    runId=requestedRunId;
    outDir=fullfile(parentDir,runId);
    suffix=1;
    while exist(outDir,'dir')
        suffix=suffix+1;
        runId=sprintf('%s_%02d',requestedRunId,suffix);
        outDir=fullfile(parentDir,runId);
    end
    [ok,msg]=mkdir(outDir);
    if ~ok
        error('diagnose_y_grid_runtime:RunDirectoryCreateFailed', ...
            'Could not create run directory %s: %s',outDir,msg);
    end
end
