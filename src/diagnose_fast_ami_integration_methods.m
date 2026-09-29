function dbundle = diagnose_fast_ami_integration_methods(source,varargin)
%DIAGNOSE_FAST_AMI_INTEGRATION_METHODS  Compare high-SNR fading quadratures.
%
% Commit I.2 diagnostic for the revised IEEE AMI pipeline.
% This function consumes a bundle produced by sim_revision_AMI_vs_SNR().
% It does NOT rerun SA and does NOT modify AMI_functions.m.
%
% For each selected validated constellation winner it compares, on the SAME
% outer y-grid/support:
%
%   (A) the existing Gauss-Hermite fading quadrature at extended orders;
%   (B) a uniform quadrature in t = ln(h), where
%           t ~ N(mu_t,sig_t^2),
%       using a finite interval mu_t +/- TSpanSigma*sig_t and trapezoidal
%       weights in t.
%
% Both approximations are scored against the independent AMI validator.
% The purpose is to decide whether the high-SNR production objective should
% continue to use higher-order GH or migrate to a resolved log-h grid.
%
% SOURCE may be either an in-memory revised_AMI_vs_SNR bundle or the path to
% revisedAMIvsSNR.mat.
%
% Name-value options:
%   'MVec'              M values ([] -> all source values)
%   'SNRVec'            SNR values ([] -> highest source SNR only)
%   'TurbulenceVec'     sigma_X^2 values ([] -> all source values)
%   'ReplicateVec'      replicate indices ([] -> all source replicates)
%   'GHVec'             extended GH orders (default [800 1200 1600 2400])
%   'DtVec'             uniform t-grid steps (default [0.02 0.01 0.005])
%   'TSpanSigma'        half-width of t interval in sigmas (default 8)
%   'YBlockSize'        y columns per log-h work block (default 512)
%   'UseParallelCases'  parallelize source-winner cases (default true)
%   'MakePlots'         create diagnostic figures (default true)
%   'SaveFigures'       save PNG/FIG pairs (default true)
%   'RunLabel'          optional label (default '')
%   'ResultsRoot'       output root (default repository/results/ieee_revision)
%
% Recommended use after Commit I.1:
%   d2 = diagnose_fast_ami_integration_methods(iPilot, ...
%       'MVec',[8 32], 'SNRVec',30, 'TurbulenceVec',0.1, ...
%       'ReplicateVec',[1 2], ...
%       'GHVec',[800 1200 1600 2400], ...
%       'DtVec',[0.02 0.01 0.005], 'TSpanSigma',8, ...
%       'UseParallelCases',true, 'RunLabel','I2_highSNR30');
%
% NOTE: GH orders above ~1000 are diagnostic only with the current dense
% Golub-Welsch implementation in AMI_functions.gaussHermite().  They may be
% computationally expensive; this experiment measures that cost explicitly.

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p,'source',@(v) isstruct(v) || ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'MVec',[],@(v) isempty(v) || numeric_vector(v));
    addParameter(p,'SNRVec',[],@(v) isempty(v) || numeric_vector(v));
    addParameter(p,'TurbulenceVec',[],@(v) isempty(v) || (numeric_vector(v) && all(v>=0)));
    addParameter(p,'ReplicateVec',[],@(v) isempty(v) || (numeric_vector(v) && all(v>=1) && all(mod(v,1)==0)));
    addParameter(p,'GHVec',[800 1200 1600 2400],@(v) numeric_vector(v) && all(v>=2) && all(mod(v,1)==0));
    addParameter(p,'DtVec',[0.02 0.01 0.005],@(v) numeric_vector(v) && all(v>0));
    addParameter(p,'TSpanSigma',8,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p,'YBlockSize',512,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=32 && mod(v,1)==0);
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
    if isempty(o.SNRVec), snrVec=max(srcSNR); else, snrVec=unique(double(o.SNRVec(:).'),'stable'); end
    if isempty(o.TurbulenceVec), sigVec=srcSig; else, sigVec=unique(double(o.TurbulenceVec(:).'),'stable'); end
    if isempty(o.ReplicateVec), repVec=srcRep; else, repVec=unique(double(o.ReplicateVec(:).'),'stable'); end
    ghVec=unique(double(o.GHVec(:).'),'stable');
    dtVec=unique(double(o.DtVec(:).'),'stable');

    idxM=map_axis(MVec,srcM,'M');
    idxSNR=map_axis(snrVec,srcSNR,'SNR');
    idxSig=map_axis(sigVec,srcSig,'sigma_X_sq');
    idxRep=map_axis(repVec,srcRep,'replicate');

    nM=numel(MVec); nS=numel(snrVec); nG=numel(sigVec); nR=numel(repVec);
    [iMGrid,iGGrid,iSGrid,iRGrid]=ndgrid(1:nM,1:nG,1:nS,1:nR);
    taskM=iMGrid(:); taskG=iGGrid(:); taskS=iSGrid(:); taskR=iRGrid(:);
    nCases=numel(taskM);

    resultsRoot=char(o.ResultsRoot);
    parentDir=fullfile(resultsRoot,'fast_ami_integration_diagnostic');
    labelToken=sanitize_label(o.RunLabel);
    ghToken=sprintf('GH%s',join_ints(ghVec));
    dtToken=sprintf('DT%s',join_decimals(dtVec));
    utcToken=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requestedRunId=sprintf('%s_%s_%s_Tspan%s_%s', ...
        labelToken,ghToken,dtToken,num_token(o.TSpanSigma,2),utcToken);
    [outDir,runId]=make_unique_run_directory(parentDir,requestedRunId);
    caseDir=fullfile(outDir,'cases'); mkdir(caseDir);
    figDir=fullfile(outDir,'figures');
    if o.MakePlots && o.SaveFigures, mkdir(figDir); end

    runMeta=fso_result_utils.run_metadata(mfilename);
    runMeta.runId=runId;
    runMeta.runLabel=char(string(o.RunLabel));
    runMeta.outputDirectory=outDir;

    fprintf('\n============================================================\n');
    fprintf('Commit I.2: high-SNR fading-integration diagnostic\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('Source: %s\n',sourceDescription);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | reps=%s\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),mat2str(repVec));
    fprintf('Extended GH orders=%s\n',mat2str(ghVec));
    fprintf('Uniform log-h dt=%s | t span=mu +/- %.2f sigma\n',mat2str(dtVec),o.TSpanSigma);
    fprintf('Cases=%d | total method evaluations=%d\n', ...
        nCases,nCases*(numel(ghVec)+numel(dtVec)));
    fprintf('Same outer y-grid used by both methods.\n');
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
            [caseResults{k},caseFiles{k}]=run_case( ...
                taskM(k),taskG(k),taskS(k),taskR(k), ...
                MVec,sigVec,snrVec,repVec,idxM,idxSig,idxSNR,idxRep, ...
                ghVec,dtVec,o.TSpanSigma,o.YBlockSize,src,runMeta,caseDir);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k=1:nCases
            [caseResults{k},caseFiles{k}]=run_case( ...
                taskM(k),taskG(k),taskS(k),taskR(k), ...
                MVec,sigVec,snrVec,repVec,idxM,idxSig,idxSNR,idxRep, ...
                ghVec,dtVec,o.TSpanSigma,o.YBlockSize,src,runMeta,caseDir);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    caseGrid=reshape(caseResults,nM,nG,nS,nR);
    fileGrid=reshape(caseFiles,nM,nG,nS,nR);
    summary=build_summary(caseGrid,ghVec,dtVec);

    dbundle=struct();
    dbundle.meta=runMeta;
    dbundle.experiment='fast_ami_integration_diagnostic';
    dbundle.run=struct('id',runId,'label',char(string(o.RunLabel)), ...
        'outputDirectory',outDir);
    dbundle.source=struct('description',sourceDescription,'runId',source_run_id(src));
    dbundle.axes=struct('M',MVec,'SNR_dB',snrVec,'sigma_X_sq',sigVec, ...
        'replicate',repVec,'GH',ghVec,'dt',dtVec);
    dbundle.settings=struct('TSpanSigma',double(o.TSpanSigma), ...
        'YBlockSize',double(o.YBlockSize),'sameOuterYGrid',true, ...
        'outerParallelism',doPar);
    dbundle.caseFiles=fileGrid;
    dbundle.summary=summary;
    dbundle.figureFiles={};

    if o.MakePlots
        dbundle.figureFiles=make_plots(dbundle,caseGrid,figDir,logical(o.SaveFigures));
    end

    bundlePath=fullfile(outDir,'fastAMIIntegrationDiagnostic.mat');
    save(bundlePath,'dbundle','-v7.3');
    fprintf('\nSaved integration diagnostic bundle: %s\n',bundlePath);
    print_summary(dbundle,caseGrid);

    function progress_tick()
        doneCount=doneCount+1;
        fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function [c,path]=run_case(iM,iG,iS,iR,MVec,sigVec,snrVec,repVec, ...
        idxM,idxSig,idxSNR,idxRep,ghVec,dtVec,tSpanSigma,yBlockSize,src,runMeta,caseDir)

    M=MVec(iM); sig=sigVec(iG); snr=snrVec(iS); rep=repVec(iR);
    sM=idxM(iM); sG=idxSig(iG); sS=idxSNR(iS); sR=idxRep(iR);

    x=real(src.summary.xWinner{sM,sG,sS,sR}(:));
    storedVal=src.summary.raw.amiGSValidated(sM,sG,sS,sR);
    storedFast=src.summary.raw.amiGSFast(sM,sG,sS,sR);

    Pavg=double(src.axes.P_avg);
    minGap=double(src.method.minGap);
    pinZero=logical(src.method.pinZero);
    xMaxBound=double(src.xMaxBounds(sM));

    cfg=build_fso_config(M,Pavg,snr,sig, ...
        'ghN_h',max(2,ghVec(1)),'xMaxBound',xMaxBound, ...
        'saMaxIter',1,'saNStarts',1,'saUseParallel',false, ...
        'minGap',minGap,'pinZero',pinZero,'seedInit',1,'logEvery',0);

    validatorAMI=cfg.AMI_Validator(x);
    validatorMismatch=validatorAMI-storedVal;
    yGrid=AMI_functions.build_noCSI_y_grid(cfg,xMaxBound);

    nGH=numel(ghVec);
    ghFastMI=nan(1,nGH); ghSignedError=nan(1,nGH); ghAbsError=nan(1,nGH); ghRuntime=nan(1,nGH);
    for i=1:nGH
        t=tic;
        v=AMI_functions.AMI_noCSI_fast_grid(x,cfg.px,cfg,ghVec(i),yGrid);
        ghRuntime(i)=toc(t);
        ghFastMI(i)=v;
        ghSignedError(i)=v-validatorAMI;
        ghAbsError(i)=abs(ghSignedError(i));
    end

    nDt=numel(dtVec);
    logtFastMI=nan(1,nDt); logtSignedError=nan(1,nDt); logtAbsError=nan(1,nDt); logtRuntime=nan(1,nDt);
    logtN=nan(1,nDt); capturedMass=nan(1,nDt); actualDtMax=nan(1,nDt);
    for i=1:nDt
        t=tic;
        [v,dg]=ami_noCSI_uniform_logh(x,cfg.px,cfg,dtVec(i),tSpanSigma,yGrid,yBlockSize);
        logtRuntime(i)=toc(t);
        logtFastMI(i)=v;
        logtSignedError(i)=v-validatorAMI;
        logtAbsError(i)=abs(logtSignedError(i));
        logtN(i)=dg.nT;
        capturedMass(i)=dg.capturedMass;
        actualDtMax(i)=dg.maxStep;
    end

    c=struct();
    c.meta=runMeta;
    c.case=struct('M',M,'SNR_dB',snr,'sigma_X_sq',sig,'replicate',rep);
    c.source=struct('x',x,'storedValidatedAMI',storedVal,'storedFastAMI',storedFast, ...
        'storedFastMinusValidated',storedFast-storedVal, ...
        'recomputedValidatorAMI',validatorAMI,'validatorMismatch',validatorMismatch, ...
        'actualXMax',max(x),'xMaxBound',xMaxBound,'xMaxToBoundRatio',max(x)/xMaxBound);
    c.outerYGrid=struct('N',numel(yGrid),'min',yGrid(1),'max',yGrid(end));
    c.GH=struct('orders',ghVec,'fastMI',ghFastMI,'signedError',ghSignedError, ...
        'absError',ghAbsError,'runtimeSeconds',ghRuntime);
    c.logH=struct('dtRequested',dtVec,'fastMI',logtFastMI,'signedError',logtSignedError, ...
        'absError',logtAbsError,'runtimeSeconds',logtRuntime,'nT',logtN, ...
        'capturedMass',capturedMass,'maxActualStep',actualDtMax, ...
        'spanSigma',tSpanSigma);

    fname=sprintf('integration_M%d_SNR%s_sig%s_rep%02d.mat', ...
        M,num_token(snr,2),num_token(sig,4),rep);
    path=fullfile(caseDir,fname);
    save(path,'c','-v7.3');
end


function [miBits,dg]=ami_noCSI_uniform_logh(x,px,params,dt,tSpanSigma,yGrid,yBlockSize)
%AMI_NOCSI_UNIFORM_LOGH  No-CSI AMI with uniform quadrature in t=ln(h).
% The outer y integration is intentionally identical to the production fast
% evaluator so this diagnostic isolates the fading quadrature.

    x=x(:).'; px=px(:).'; yGrid=yGrid(:).';
    M=numel(x); Ny=numel(yGrid);

    if params.sig_t < 1e-12
        % Degenerate no-turbulence case: the fading integral disappears.
        means=params.R*x(:);
        d2=bsxfun(@minus,yGrid,means).^2;
        logPyx=-0.5*log(2*pi*params.sigma_n_sq)-d2/(2*params.sigma_n_sq);
        miBits=mi_from_log_pyx(logPyx,px,yGrid);
        dg=struct('nT',1,'capturedMass',1,'maxStep',0, ...
            'tMin',params.mu_t,'tMax',params.mu_t);
        return;
    end

    tMin=params.mu_t-tSpanSigma*params.sig_t;
    tMax=params.mu_t+tSpanSigma*params.sig_t;
    tGrid=tMin:dt:tMax;
    if isempty(tGrid) || tGrid(end)<tMax-100*eps(max(1,abs(tMax)))
        tGrid=[tGrid tMax]; %#ok<AGROW>
    else
        tGrid(end)=tMax;
    end
    tGrid=unique(tGrid,'stable');

    qWeights=trapezoid_node_weights(tGrid);
    fT=exp(-0.5*((tGrid-params.mu_t)/params.sig_t).^2) / ...
        (sqrt(2*pi)*params.sig_t);
    q=qWeights.*fT;
    capturedMass=sum(q);

    % Zero/underflow quadrature masses do not contribute to the mixture.
    keep=q>0 & isfinite(q);
    tGrid=tGrid(keep); q=q(keep);
    logQ=log(q(:));
    hVals=exp(tGrid(:));

    inv2s2=1/(2*params.sigma_n_sq);
    logNorm=-0.5*log(2*pi*params.sigma_n_sq);
    logPyx=zeros(M,Ny);

    for j=1:M
        means=params.R*hVals*x(j);
        for b0=1:yBlockSize:Ny
            b1=min(Ny,b0+yBlockSize-1);
            yb=yGrid(b0:b1);
            d2=bsxfun(@minus,yb,means).^2;
            logTerms=bsxfun(@plus,logQ+logNorm,-inv2s2*d2);
            logPyx(j,b0:b1)=local_logsumexp(logTerms,1);
        end
    end

    miBits=mi_from_log_pyx(logPyx,px,yGrid);
    steps=diff(tGrid);
    if isempty(steps), maxStep=0; else, maxStep=max(steps); end
    dg=struct('nT',numel(tGrid),'capturedMass',capturedMass,'maxStep',maxStep, ...
        'tMin',tMin,'tMax',tMax);
end


function miBits=mi_from_log_pyx(logPyx,px,yGrid)
    logPxPyx=bsxfun(@plus,log(px(:)),logPyx);
    logPy=local_logsumexp(logPxPyx,1);
    M=size(logPyx,1);
    MI=0;
    for j=1:M
        Pyxj=exp(logPyx(j,:));
        lr=(logPyx(j,:)-logPy)/log(2);
        integrand=Pyxj.*lr;
        integrand(~isfinite(integrand))=0;
        MI=MI+px(j)*trapz(yGrid,integrand);
    end
    miBits=max(0,MI);
end


function w=trapezoid_node_weights(x)
    x=x(:).'; n=numel(x);
    if n==1, w=1; return; end
    dx=diff(x);
    w=zeros(1,n);
    w(1)=dx(1)/2;
    w(end)=dx(end)/2;
    if n>2
        w(2:end-1)=(dx(1:end-1)+dx(2:end))/2;
    end
end


function s=local_logsumexp(A,dim)
    if nargin<2, dim=1; end
    m=max(A,[],dim);
    s=m+log(sum(exp(bsxfun(@minus,A,m)),dim));
end


function S=build_summary(caseGrid,ghVec,dtVec)
    nGH=numel(ghVec); nDt=numel(dtVec);
    ghAbs=nan(nGH,numel(caseGrid)); ghSigned=ghAbs; ghRun=ghAbs;
    dtAbs=nan(nDt,numel(caseGrid)); dtSigned=dtAbs; dtRun=dtAbs; dtN=dtAbs; dtMass=dtAbs;
    valMismatch=nan(1,numel(caseGrid)); ratios=valMismatch;

    for k=1:numel(caseGrid)
        c=caseGrid{k};
        ghAbs(:,k)=c.GH.absError(:);
        ghSigned(:,k)=c.GH.signedError(:);
        ghRun(:,k)=c.GH.runtimeSeconds(:);
        dtAbs(:,k)=c.logH.absError(:);
        dtSigned(:,k)=c.logH.signedError(:);
        dtRun(:,k)=c.logH.runtimeSeconds(:);
        dtN(:,k)=c.logH.nT(:);
        dtMass(:,k)=c.logH.capturedMass(:);
        valMismatch(k)=c.source.validatorMismatch;
        ratios(k)=c.source.xMaxToBoundRatio;
    end

    S=struct();
    S.GH=struct();
    S.GH.orders=ghVec;
    S.GH.meanAbsError=mean(ghAbs,2,'omitnan').';
    S.GH.maxAbsError=max(ghAbs,[],2,'omitnan').';
    S.GH.meanSignedError=mean(ghSigned,2,'omitnan').';
    S.GH.meanRuntimeSeconds=mean(ghRun,2,'omitnan').';
    S.GH.maxRuntimeSeconds=max(ghRun,[],2,'omitnan').';

    S.logH=struct();
    S.logH.dt=dtVec;
    S.logH.meanAbsError=mean(dtAbs,2,'omitnan').';
    S.logH.maxAbsError=max(dtAbs,[],2,'omitnan').';
    S.logH.meanSignedError=mean(dtSigned,2,'omitnan').';
    S.logH.meanRuntimeSeconds=mean(dtRun,2,'omitnan').';
    S.logH.maxRuntimeSeconds=max(dtRun,[],2,'omitnan').';
    S.logH.meanNodes=mean(dtN,2,'omitnan').';
    S.logH.minCapturedMass=min(dtMass,[],2,'omitnan').';

    S.maxAbsValidatorRecomputeMismatch=max(abs(valMismatch),[],'omitnan');
    S.meanXMaxToBoundRatio=mean(ratios,'omitnan');
    S.maxXMaxToBoundRatio=max(ratios,[],'omitnan');
    S.thresholds=struct('absError1e2',1e-2,'absError1e3',1e-3,'absError1e4',1e-4);
    S.GH.firstOrderAllCasesLe1e2=first_setting(ghVec,S.GH.maxAbsError,1e-2,'ascending');
    S.GH.firstOrderAllCasesLe1e3=first_setting(ghVec,S.GH.maxAbsError,1e-3,'ascending');
    S.GH.firstOrderAllCasesLe1e4=first_setting(ghVec,S.GH.maxAbsError,1e-4,'ascending');
    S.logH.coarsestDtAllCasesLe1e2=first_setting(dtVec,S.logH.maxAbsError,1e-2,'descending');
    S.logH.coarsestDtAllCasesLe1e3=first_setting(dtVec,S.logH.maxAbsError,1e-3,'descending');
    S.logH.coarsestDtAllCasesLe1e4=first_setting(dtVec,S.logH.maxAbsError,1e-4,'descending');
end


function v=first_setting(settings,errors,tol,direction)
    settings=settings(:).'; errors=errors(:).';
    if strcmp(direction,'ascending')
        [ss,ord]=sort(settings,'ascend');
    else
        [ss,ord]=sort(settings,'descend');
    end
    ee=errors(ord);
    idx=find(ee<=tol,1,'first');
    if isempty(idx), v=NaN; else, v=ss(idx); end
end


function files=make_plots(bundle,caseGrid,figDir,saveFigures)
    files={};
    S=bundle.summary;

    f1=figure('Name','High-SNR integration accuracy','Color','w');
    tl=tiledlayout(f1,1,2,'TileSpacing','compact','Padding','compact');
    title(tl,'Fast AMI fading-integration accuracy vs independent validator');

    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    semilogy(ax,S.GH.orders,S.GH.meanAbsError,'-o','LineWidth',1.5,'DisplayName','mean');
    semilogy(ax,S.GH.orders,S.GH.maxAbsError,'-s','LineWidth',1.5,'DisplayName','worst case');
    yline(ax,1e-3,'--','10^{-3}','HandleVisibility','off');
    xlabel(ax,'Gauss-Hermite order'); ylabel(ax,'|I_{fast}-I_{val}| [bit/symbol]');
    title(ax,'Extended Gauss-Hermite'); legend(ax,'Location','best');

    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    semilogy(ax,S.logH.meanNodes,S.logH.meanAbsError,'-o','LineWidth',1.5,'DisplayName','mean');
    semilogy(ax,S.logH.meanNodes,S.logH.maxAbsError,'-s','LineWidth',1.5,'DisplayName','worst case');
    yline(ax,1e-3,'--','10^{-3}','HandleVisibility','off');
    xlabel(ax,'Uniform log-h nodes'); ylabel(ax,'|I_{fast}-I_{val}| [bit/symbol]');
    title(ax,'Uniform t=ln(h) quadrature'); legend(ax,'Location','best');

    f2=figure('Name','High-SNR integration runtime','Color','w');
    tl=tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
    title(tl,'Evaluation cost of candidate fading quadratures');
    ax=nexttile(tl); grid(ax,'on'); box(ax,'on');
    plot(ax,S.GH.orders,S.GH.meanRuntimeSeconds,'-o','LineWidth',1.5);
    xlabel(ax,'Gauss-Hermite order'); ylabel(ax,'Mean AMI evaluation time [s]'); title(ax,'Extended GH');
    ax=nexttile(tl); grid(ax,'on'); box(ax,'on');
    plot(ax,S.logH.meanNodes,S.logH.meanRuntimeSeconds,'-o','LineWidth',1.5);
    xlabel(ax,'Uniform log-h nodes'); ylabel(ax,'Mean AMI evaluation time [s]'); title(ax,'Uniform log-h');

    f3=figure('Name','Per-case integration error','Color','w');
    nCases=numel(caseGrid); nCols=min(2,nCases); nRows=ceil(nCases/nCols);
    tl=tiledlayout(f3,nRows,nCols,'TileSpacing','compact','Padding','compact');
    title(tl,'Per-constellation high-SNR integration error');
    for k=1:nCases
        c=caseGrid{k};
        ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        semilogy(ax,c.GH.orders,c.GH.absError,'-o','LineWidth',1.3,'DisplayName','GH');
        % Plot log-h against an auxiliary topological index to keep both on one panel.
        x2=linspace(min(c.GH.orders),max(c.GH.orders),numel(c.logH.absError));
        semilogy(ax,x2,c.logH.absError,'-s','LineWidth',1.3,'DisplayName','log-h grid');
        yline(ax,1e-3,'--','HandleVisibility','off');
        xlabel(ax,'GH order / log-h diagnostic index'); ylabel(ax,'absolute AMI error');
        title(ax,sprintf('M=%d, SNR=%.1f, rep=%d',c.case.M,c.case.SNR_dB,c.case.replicate));
        legend(ax,'Location','best');
    end

    if saveFigures
        if ~exist(figDir,'dir'), mkdir(figDir); end
        files=[files save_figure_pair(f1,figDir,'integration_accuracy')]; %#ok<AGROW>
        files=[files save_figure_pair(f2,figDir,'integration_runtime')]; %#ok<AGROW>
        files=[files save_figure_pair(f3,figDir,'integration_per_case')]; %#ok<AGROW>
    end
end


function print_summary(bundle,caseGrid)
    S=bundle.summary;
    fprintf('\nFAST AMI INTEGRATION-METHOD DIAGNOSTIC SUMMARY\n');
    fprintf('Extended GH orders: %s\n',mat2str(S.GH.orders));
    fprintf('Uniform log-h dt: %s\n',mat2str(S.logH.dt));
    fprintf('Max validator recompute mismatch: %.3e bit/symbol\n',S.maxAbsValidatorRecomputeMismatch);

    for k=1:numel(caseGrid)
        c=caseGrid{k};
        fprintf('\nM=%d | SNR=%.1f dB | sigma_X^2=%.4f | rep=%d\n', ...
            c.case.M,c.case.SNR_dB,c.case.sigma_X_sq,c.case.replicate);
        fprintf('  validator=%.9f | stored fast-val=%+.3e\n', ...
            c.source.recomputedValidatorAMI,c.source.storedFastMinusValidated);
        fprintf('  outer y-grid N=%d, range=[%.4f, %.4f]\n', ...
            c.outerYGrid.N,c.outerYGrid.min,c.outerYGrid.max);
        fprintf('  GH: order | abs error | runtime[s]\n');
        for i=1:numel(c.GH.orders)
            fprintf('      %5d | %.6g | %.4f\n',c.GH.orders(i),c.GH.absError(i),c.GH.runtimeSeconds(i));
        end
        fprintf('  log-h: dt | nodes | captured mass | abs error | runtime[s]\n');
        for i=1:numel(c.logH.dtRequested)
            fprintf('      %.6g | %5d | %.12f | %.6g | %.4f\n', ...
                c.logH.dtRequested(i),c.logH.nT(i),c.logH.capturedMass(i), ...
                c.logH.absError(i),c.logH.runtimeSeconds(i));
        end
    end

    fprintf('\nAggregate GH results:\n');
    fprintf('orders = %s\n',mat2str(S.GH.orders));
    fprintf('mean |error| = %s\n',mat2str(S.GH.meanAbsError,6));
    fprintf('max  |error| = %s\n',mat2str(S.GH.maxAbsError,6));
    fprintf('mean runtime = %s s\n',mat2str(S.GH.meanRuntimeSeconds,5));

    fprintf('\nAggregate uniform log-h results:\n');
    fprintf('dt = %s\n',mat2str(S.logH.dt));
    fprintf('mean nodes = %s\n',mat2str(S.logH.meanNodes,6));
    fprintf('mean |error| = %s\n',mat2str(S.logH.meanAbsError,6));
    fprintf('max  |error| = %s\n',mat2str(S.logH.maxAbsError,6));
    fprintf('mean runtime = %s s\n',mat2str(S.logH.meanRuntimeSeconds,5));
    fprintf('minimum captured probability mass = %s\n',mat2str(S.logH.minCapturedMass,12));

    fprintf('\nAll-case thresholds:\n');
    fprintf('  GH order for |error|<=1e-2: %s\n',fmt_setting(S.GH.firstOrderAllCasesLe1e2));
    fprintf('  GH order for |error|<=1e-3: %s\n',fmt_setting(S.GH.firstOrderAllCasesLe1e3));
    fprintf('  GH order for |error|<=1e-4: %s\n',fmt_setting(S.GH.firstOrderAllCasesLe1e4));
    fprintf('  coarsest log-h dt for |error|<=1e-2: %s\n',fmt_setting(S.logH.coarsestDtAllCasesLe1e2));
    fprintf('  coarsest log-h dt for |error|<=1e-3: %s\n',fmt_setting(S.logH.coarsestDtAllCasesLe1e3));
    fprintf('  coarsest log-h dt for |error|<=1e-4: %s\n',fmt_setting(S.logH.coarsestDtAllCasesLe1e4));
end


function s=fmt_setting(v)
    if isnan(v), s='not reached'; else, s=sprintf('%.6g',v); end
end


function files=save_figure_pair(f,figDir,name)
    png=fullfile(figDir,[name '.png']);
    fig=fullfile(figDir,[name '.fig']);
    exportgraphics(f,png,'Resolution',220);
    savefig(f,fig);
    files={png,fig};
end


function [src,desc]=load_source(source)
    if isstruct(source)
        src=source; desc='<workspace bundle>'; return;
    end
    path=char(string(source));
    if ~exist(path,'file')
        error('diagnose_fast_ami_integration_methods:MissingSource','Source file not found: %s',path);
    end
    D=load(path);
    if isfield(D,'bundle')
        src=D.bundle;
    elseif isfield(D,'iPilot')
        src=D.iPilot;
    else
        names=fieldnames(D); src=[];
        for k=1:numel(names)
            if isstruct(D.(names{k})) && isfield(D.(names{k}),'experiment')
                src=D.(names{k}); break;
            end
        end
        if isempty(src)
            error('diagnose_fast_ami_integration_methods:BadSourceFile','No compatible bundle found in %s.',path);
        end
    end
    desc=path;
end


function validate_source(src)
    if ~isfield(src,'experiment') || ~strcmp(src.experiment,'revised_AMI_vs_SNR')
        error('diagnose_fast_ami_integration_methods:BadSource', ...
            'Source must be a revised_AMI_vs_SNR bundle from sim_revision_AMI_vs_SNR.');
    end
    req={'axes','method','summary','xMaxBounds'};
    for k=1:numel(req)
        if ~isfield(src,req{k}), error('diagnose_fast_ami_integration_methods:BadSource','Missing source field %s.',req{k}); end
    end
    if ~isfield(src.summary,'xWinner') || ~isfield(src.summary,'raw') || ...
            ~isfield(src.summary.raw,'amiGSValidated') || ~isfield(src.summary.raw,'amiGSFast')
        error('diagnose_fast_ami_integration_methods:BadSource','Source summary lacks winner or AMI fields.');
    end
end


function idx=map_axis(requested,available,name)
    idx=nan(size(requested));
    for k=1:numel(requested)
        tol=1e-10*max(1,abs(requested(k)));
        j=find(abs(available-requested(k))<=tol,1);
        if isempty(j)
            error('diagnose_fast_ami_integration_methods:MissingAxis', ...
                'Requested %s value %.12g is not present in source.',name,requested(k));
        end
        idx(k)=j;
    end
end


function tf=numeric_vector(v)
    tf=isnumeric(v) && isvector(v) && ~isempty(v) && all(isfinite(v)) && all(isreal(v));
end


function id=source_run_id(src)
    if isfield(src,'run') && isfield(src.run,'id'), id=src.run.id; else, id=''; end
end


function token=sanitize_label(label)
    token=char(string(label)); token=strtrim(token);
    if isempty(token), token='run'; end
    token=regexprep(token,'[^A-Za-z0-9_-]+','_');
    token=regexprep(token,'_+','_'); token=regexprep(token,'^_+|_+$','');
    if isempty(token), token='run'; end
end


function s=join_ints(v)
    c=arrayfun(@(x) sprintf('%d',x),v,'UniformOutput',false);
    s=strjoin(c,'-');
end


function s=join_decimals(v)
    c=arrayfun(@(x) num_token(x,4),v,'UniformOutput',false);
    s=strjoin(c,'-');
end


function token=num_token(v,nDec)
    token=sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token=strrep(token,'-','m'); token=strrep(token,'+','p'); token=strrep(token,'.','p');
end


function [outDir,runId]=make_unique_run_directory(parentDir,requestedRunId)
    if ~exist(parentDir,'dir'), mkdir(parentDir); end
    runId=requestedRunId; outDir=fullfile(parentDir,runId); suffix=1;
    while exist(outDir,'dir')
        runId=sprintf('%s_%02d',requestedRunId,suffix);
        outDir=fullfile(parentDir,runId); suffix=suffix+1;
    end
    mkdir(outDir);
end
