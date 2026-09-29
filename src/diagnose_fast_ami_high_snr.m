function dbundle = diagnose_fast_ami_high_snr(source,varargin)
%DIAGNOSE_FAST_AMI_HIGH_SNR  Separate GH-order and y-grid errors at high SNR.
%
% Commit I.1 diagnostic for the revised AMI-vs-SNR production pipeline.
% This function consumes a saved/in-memory bundle produced by
% sim_revision_AMI_vs_SNR().  It does NOT rerun SA and does NOT modify the
% source bundle.  For each selected validated winner it recomputes:
%
%   * the independent validator AMI (reference);
%   * fast AMI for several Gauss-Hermite orders;
%   * fast AMI on the original y-grid and uniformly subdivided versions of
%     that SAME grid/support.
%
% Therefore:
%   - changing GH order at grid factor 1 isolates GH quadrature sensitivity;
%   - changing GridRefineVec at fixed GH isolates outer-y-grid sensitivity;
%   - the full GH x grid-factor matrix reveals interactions between them.
%
% SOURCE may be either an in-memory revised_AMI_vs_SNR bundle or the path to
% revisedAMIvsSNR.mat.
%
% Name-value options:
%   'MVec'              M values ([] -> all source values)
%   'SNRVec'            SNR values ([] -> highest source SNR only)
%   'TurbulenceVec'     sigma_X^2 values ([] -> all source values)
%   'ReplicateVec'      replicate indices ([] -> all source replicates)
%   'GHVec'             GH orders (default [200 400 600 800])
%   'GridRefineVec'     integer interval subdivision factors (default [1 2])
%   'UseParallelCases'  parallelize source-winner cases (default true)
%   'MakePlots'         create diagnostic figure (default true)
%   'SaveFigures'       save PNG/FIG pair (default true)
%   'RunLabel'          optional label (default '')
%   'ResultsRoot'       output root (default repository/results/ieee_revision)
%
% Recommended first high-SNR diagnostic:
%   d = diagnose_fast_ami_high_snr(iPilot, ...
%       'MVec',[8 32], 'SNRVec',30, 'TurbulenceVec',0.1, ...
%       'ReplicateVec',[1 2], 'GHVec',[200 400 600 800], ...
%       'GridRefineVec',[1 2], 'UseParallelCases',true, ...
%       'RunLabel','I1_highSNR30');
%
% Interpretation:
%   If error falls strongly with GH at GridRefine=1, GH resolution is the
%   dominant problem.  If it barely changes with GH but falls when the grid
%   is refined, the y-grid is dominant.  If both are required, the two
%   discretizations interact.

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p,'source',@(v) isstruct(v) || ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'MVec',[],@(v) isempty(v) || numeric_vector(v));
    addParameter(p,'SNRVec',[],@(v) isempty(v) || numeric_vector(v));
    addParameter(p,'TurbulenceVec',[],@(v) isempty(v) || (numeric_vector(v) && all(v>=0)));
    addParameter(p,'ReplicateVec',[],@(v) isempty(v) || (numeric_vector(v) && all(v>=1) && all(mod(v,1)==0)));
    addParameter(p,'GHVec',[200 400 600 800],@(v) numeric_vector(v) && all(v>=2) && all(mod(v,1)==0));
    addParameter(p,'GridRefineVec',[1 2],@(v) numeric_vector(v) && all(v>=1) && all(mod(v,1)==0));
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

    if isempty(o.MVec), MVec = srcM; else, MVec = unique(double(o.MVec(:).'),'stable'); end
    if isempty(o.SNRVec), snrVec = max(srcSNR); else, snrVec = unique(double(o.SNRVec(:).'),'stable'); end
    if isempty(o.TurbulenceVec), sigVec = srcSig; else, sigVec = unique(double(o.TurbulenceVec(:).'),'stable'); end
    if isempty(o.ReplicateVec), repVec = srcRep; else, repVec = unique(double(o.ReplicateVec(:).'),'stable'); end
    ghVec = unique(double(o.GHVec(:).'),'stable');
    refineVec = unique(double(o.GridRefineVec(:).'),'stable');

    idxM = map_axis(MVec,srcM,'M');
    idxSNR = map_axis(snrVec,srcSNR,'SNR');
    idxSig = map_axis(sigVec,srcSig,'sigma_X_sq');
    idxRep = map_axis(repVec,srcRep,'replicate');

    nM=numel(MVec); nS=numel(snrVec); nG=numel(sigVec); nR=numel(repVec);
    [iMGrid,iGGrid,iSGrid,iRGrid] = ndgrid(1:nM,1:nG,1:nS,1:nR);
    taskM=iMGrid(:); taskG=iGGrid(:); taskS=iSGrid(:); taskR=iRGrid(:);
    nCases=numel(taskM);

    resultsRoot = char(o.ResultsRoot);
    parentDir = fullfile(resultsRoot,'fast_ami_high_snr_diagnostic');
    labelToken = sanitize_label(o.RunLabel);
    ghToken = sprintf('GH%s',join_ints(ghVec));
    gridToken = sprintf('Yx%s',join_ints(refineVec));
    utcToken = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requestedRunId = sprintf('%s_%s_%s_%s',labelToken,ghToken,gridToken,utcToken);
    [outDir,runId] = make_unique_run_directory(parentDir,requestedRunId);
    caseDir = fullfile(outDir,'cases'); mkdir(caseDir);
    figDir = fullfile(outDir,'figures');
    if o.MakePlots && o.SaveFigures, mkdir(figDir); end

    runMeta = fso_result_utils.run_metadata(mfilename);
    runMeta.runId = runId;
    runMeta.runLabel = char(string(o.RunLabel));
    runMeta.outputDirectory = outDir;

    fprintf('\n============================================================\n');
    fprintf('Commit I.1: high-SNR fast-AMI diagnostic\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('Source: %s\n',sourceDescription);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | reps=%s\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),mat2str(repVec));
    fprintf('GH orders=%s | y-grid refine factors=%s\n',mat2str(ghVec),mat2str(refineVec));
    fprintf('Cases=%d | fast evaluations=%d\n',nCases,nCases*numel(ghVec)*numel(refineVec));
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    caseResults = cell(nCases,1);
    caseFiles = cell(nCases,1);

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
        if isempty(pool), pool=parpool('Processes'); end
        fprintf('Outer case parallelism: %d process workers\n',pool.NumWorkers);
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~) progress_tick());
        parfor k=1:nCases
            [caseResults{k},caseFiles{k}] = run_case( ...
                taskM(k),taskG(k),taskS(k),taskR(k), ...
                MVec,sigVec,snrVec,repVec,idxM,idxSig,idxSNR,idxRep, ...
                ghVec,refineVec,src,runMeta,caseDir);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k=1:nCases
            [caseResults{k},caseFiles{k}] = run_case( ...
                taskM(k),taskG(k),taskS(k),taskR(k), ...
                MVec,sigVec,snrVec,repVec,idxM,idxSig,idxSNR,idxRep, ...
                ghVec,refineVec,src,runMeta,caseDir);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    caseGrid = reshape(caseResults,nM,nG,nS,nR);
    fileGrid = reshape(caseFiles,nM,nG,nS,nR);
    summary = build_summary(caseGrid,ghVec,refineVec);

    dbundle = struct();
    dbundle.meta = runMeta;
    dbundle.experiment = 'fast_ami_high_snr_diagnostic';
    dbundle.run = struct('id',runId,'label',char(string(o.RunLabel)), ...
        'outputDirectory',outDir);
    dbundle.source = struct('description',sourceDescription, ...
        'runId',source_run_id(src));
    dbundle.axes = struct('M',MVec,'SNR_dB',snrVec,'sigma_X_sq',sigVec, ...
        'replicate',repVec,'GH',ghVec,'gridRefineFactor',refineVec);
    dbundle.caseFiles = fileGrid;
    dbundle.summary = summary;
    dbundle.figureFiles = {};

    if o.MakePlots
        dbundle.figureFiles = make_plots(dbundle,caseGrid,figDir,logical(o.SaveFigures));
    end

    bundlePath = fullfile(outDir,'fastAMIHighSNRDiagnostic.mat');
    save(bundlePath,'dbundle','-v7.3');
    fprintf('\nSaved diagnostic bundle: %s\n',bundlePath);
    print_summary(dbundle,caseGrid);

    function progress_tick()
        doneCount=doneCount+1;
        fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function [c,path] = run_case(iM,iG,iS,iR,MVec,sigVec,snrVec,repVec, ...
        idxM,idxSig,idxSNR,idxRep,ghVec,refineVec,src,runMeta,caseDir)

    M=MVec(iM); sig=sigVec(iG); snr=snrVec(iS); rep=repVec(iR);
    sM=idxM(iM); sG=idxSig(iG); sS=idxSNR(iS); sR=idxRep(iR);

    x = real(src.summary.xWinner{sM,sG,sS,sR}(:));
    storedVal = src.summary.raw.amiGSValidated(sM,sG,sS,sR);
    storedFast = src.summary.raw.amiGSFast(sM,sG,sS,sR);

    Pavg = double(src.axes.P_avg);
    minGap = double(src.method.minGap);
    pinZero = logical(src.method.pinZero);
    xMaxBound = double(src.xMaxBounds(sM));

    cfg = build_fso_config(M,Pavg,snr,sig, ...
        'ghN_h',ghVec(1),'xMaxBound',xMaxBound, ...
        'saMaxIter',1,'saNStarts',1,'saUseParallel',false, ...
        'minGap',minGap,'pinZero',pinZero,'seedInit',1,'logEvery',0);

    validatorAMI = cfg.AMI_Validator(x);
    validatorMismatch = validatorAMI-storedVal;

    baseGrid = AMI_functions.build_noCSI_y_grid(cfg,xMaxBound);
    nGH=numel(ghVec); nF=numel(refineVec);
    fastMI=nan(nGH,nF); signedError=nan(nGH,nF); absError=nan(nGH,nF);
    runtimeSeconds=nan(nGH,nF); gridN=nan(1,nF);

    grids=cell(1,nF);
    for j=1:nF
        grids{j}=refine_grid(baseGrid,refineVec(j));
        gridN(j)=numel(grids{j});
    end

    for j=1:nF
        yg=grids{j};
        for i=1:nGH
            t=tic;
            v=AMI_functions.AMI_noCSI_fast_grid(x,cfg.px,cfg,ghVec(i),yg);
            runtimeSeconds(i,j)=toc(t);
            fastMI(i,j)=v;
            signedError(i,j)=v-validatorAMI;
            absError(i,j)=abs(signedError(i,j));
        end
    end

    c=struct();
    c.meta=runMeta;
    c.case=struct('M',M,'SNR_dB',snr,'sigma_X_sq',sig,'replicate',rep);
    c.source=struct('x',x,'storedValidatedAMI',storedVal,'storedFastAMI',storedFast, ...
        'storedFastMinusValidated',storedFast-storedVal, ...
        'recomputedValidatorAMI',validatorAMI,'validatorMismatch',validatorMismatch, ...
        'actualXMax',max(x),'xMaxBound',xMaxBound, ...
        'xMaxToBoundRatio',max(x)/xMaxBound);
    c.grid=struct('baseN',numel(baseGrid),'baseMin',baseGrid(1),'baseMax',baseGrid(end), ...
        'refineFactors',refineVec,'N',gridN);
    c.GH=ghVec;
    c.fastMI=fastMI;
    c.signedError=signedError;
    c.absError=absError;
    c.runtimeSeconds=runtimeSeconds;

    fname=sprintf('diag_M%d_SNR%s_sig%s_rep%02d.mat', ...
        M,num_token(snr,2),num_token(sig,4),rep);
    path=fullfile(caseDir,fname);
    save(path,'c','-v7.3');
end


function S = build_summary(caseGrid,ghVec,refineVec)
    nGH=numel(ghVec); nF=numel(refineVec);
    allAbs=[]; allSigned=[]; allRuntime=[]; storedGap=[]; valMismatch=[]; ratios=[];
    for k=1:numel(caseGrid)
        c=caseGrid{k};
        allAbs=cat(3,allAbs,c.absError);
        allSigned=cat(3,allSigned,c.signedError);
        allRuntime=cat(3,allRuntime,c.runtimeSeconds);
        storedGap(end+1)=c.source.storedFastMinusValidated; %#ok<AGROW>
        valMismatch(end+1)=c.source.validatorMismatch; %#ok<AGROW>
        ratios(end+1)=c.source.xMaxToBoundRatio; %#ok<AGROW>
    end
    S=struct();
    S.GH=ghVec; S.gridRefineFactor=refineVec;
    S.meanAbsError=mean(allAbs,3,'omitnan');
    S.maxAbsError=max(allAbs,[],3,'omitnan');
    S.meanSignedError=mean(allSigned,3,'omitnan');
    S.meanRuntimeSeconds=mean(allRuntime,3,'omitnan');
    S.maxRuntimeSeconds=max(allRuntime,[],3,'omitnan');
    S.maxAbsStoredFastValidationGap=max(abs(storedGap),[],'omitnan');
    S.maxAbsValidatorRecomputeMismatch=max(abs(valMismatch),[],'omitnan');
    S.meanXMaxToBoundRatio=mean(ratios,'omitnan');
    S.maxXMaxToBoundRatio=max(ratios,[],'omitnan');
end


function files = make_plots(bundle,caseGrid,figDir,saveFigures)
    gh=bundle.axes.GH; rf=bundle.axes.gridRefineFactor; files={};
    f=figure('Name','High-SNR fast AMI diagnostic','Color','w');
    nCases=numel(caseGrid); nCols=min(2,nCases); nRows=ceil(nCases/nCols);
    tl=tiledlayout(f,nRows,nCols,'TileSpacing','compact','Padding','compact');
    title(tl,'Fast AMI error vs independent validator');
    for k=1:nCases
        c=caseGrid{k}; ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for j=1:numel(rf)
            semilogy(ax,gh,max(c.absError(:,j),realmin),'-o','LineWidth',1.3, ...
                'DisplayName',sprintf('y-grid x%d',rf(j)));
        end
        xlabel(ax,'Gauss-Hermite order'); ylabel(ax,'|I_{fast}-I_{val}| [bit/symbol]');
        title(ax,sprintf('M=%d, SNR=%.1f dB, rep=%d', ...
            c.case.M,c.case.SNR_dB,c.case.replicate));
        legend(ax,'Location','best');
    end
    if saveFigures
        files=save_figure_pair(f,figDir,'fast_ami_high_snr_error');
    end
end


function print_summary(bundle,caseGrid)
    fprintf('\nFAST AMI HIGH-SNR DIAGNOSTIC SUMMARY\n');
    fprintf('GH orders: %s\n',mat2str(bundle.axes.GH));
    fprintf('Grid refine factors: %s\n',mat2str(bundle.axes.gridRefineFactor));
    fprintf('Max validator recompute mismatch: %.3e bit/symbol\n', ...
        bundle.summary.maxAbsValidatorRecomputeMismatch);

    for k=1:numel(caseGrid)
        c=caseGrid{k};
        fprintf('\nM=%d | SNR=%.1f dB | sigma_X^2=%.4f | rep=%d\n', ...
            c.case.M,c.case.SNR_dB,c.case.sigma_X_sq,c.case.replicate);
        fprintf('  validator=%.9f | stored fast-val=%+.3e\n', ...
            c.source.recomputedValidatorAMI,c.source.storedFastMinusValidated);
        fprintf('  xMax=%.4f | feasible xMax bound=%.4f | ratio=%.3f\n', ...
            c.source.actualXMax,c.source.xMaxBound,c.source.xMaxToBoundRatio);
        fprintf('  base y-grid N=%d, range=[%.4f, %.4f]\n', ...
            c.grid.baseN,c.grid.baseMin,c.grid.baseMax);
        fprintf('  abs error rows=GH, cols=grid factors %s\n',mat2str(c.grid.refineFactors));
        for i=1:numel(c.GH)
            fprintf('    GH=%4d : %s\n',c.GH(i),mat2str(c.absError(i,:),6));
        end
    end

    fprintf('\nWorst-case |fast-validator| across selected cases:\n');
    fprintf('rows GH=%s, cols grid factors=%s\n', ...
        mat2str(bundle.axes.GH),mat2str(bundle.axes.gridRefineFactor));
    disp(bundle.summary.maxAbsError);
    fprintf('Mean |fast-validator| across selected cases:\n');
    disp(bundle.summary.meanAbsError);
    fprintf('Mean evaluation runtime [s]:\n');
    disp(bundle.summary.meanRuntimeSeconds);
end


function y2 = refine_grid(y,factor)
    y=double(y(:).');
    factor=double(factor);
    if factor==1
        y2=y;
        return;
    end
    n=numel(y);
    u=1:n;
    uq=linspace(1,n,(n-1)*factor+1);
    y2=interp1(u,y,uq,'linear');
    y2=unique(y2);
end


function [src,description] = load_source(source)
    if isstruct(source)
        src=source; description='<workspace bundle>'; return;
    end
    path=char(string(source));
    if ~isfile(path), error('diagnose_fast_ami_high_snr:SourceNotFound','Source file not found: %s',path); end
    D=load(path);
    if isfield(D,'bundle')
        src=D.bundle;
    elseif isfield(D,'dbundle')
        error('diagnose_fast_ami_high_snr:WrongSource','Supplied file is already a diagnostic bundle.');
    else
        fn=fieldnames(D);
        if numel(fn)==1 && isstruct(D.(fn{1})), src=D.(fn{1});
        else, error('diagnose_fast_ami_high_snr:BundleMissing','Could not identify source bundle.'); end
    end
    description=path;
end


function validate_source(src)
    if ~isfield(src,'experiment') || ~strcmp(src.experiment,'revised_AMI_vs_SNR')
        error('diagnose_fast_ami_high_snr:BadSource','Source must be a revised_AMI_vs_SNR bundle.');
    end
    required={'axes','method','summary','xMaxBounds'};
    for k=1:numel(required)
        if ~isfield(src,required{k}), error('diagnose_fast_ami_high_snr:BadSource','Missing source field %s.',required{k}); end
    end
    if ~isfield(src.summary,'xWinner') || ~isfield(src.summary,'raw')
        error('diagnose_fast_ami_high_snr:BadSource','Source summary lacks xWinner/raw fields.');
    end
end


function idx = map_axis(req,avail,name)
    idx=zeros(size(req)); scale=max(1,max(abs(avail))); tol=1e-10*scale;
    for k=1:numel(req)
        j=find(abs(avail-req(k))<=tol,1,'first');
        if isempty(j), error('diagnose_fast_ami_high_snr:AxisValueMissing','Requested %s %.12g is absent.',name,req(k)); end
        idx(k)=j;
    end
end

function id=source_run_id(src)
    id=''; if isfield(src,'run') && isfield(src.run,'id'), id=src.run.id; end
end

function tf=numeric_vector(v)
    tf=isnumeric(v) && isvector(v) && ~isempty(v) && all(isfinite(v)) && all(isreal(v));
end

function token=join_ints(v)
    parts=arrayfun(@(x) sprintf('%d',x),v,'UniformOutput',false); token=strjoin(parts,'-');
end

function token=sanitize_label(v)
    token=strtrim(char(string(v))); if isempty(token), token='run'; return; end
    token=regexprep(token,'[^A-Za-z0-9._-]+','_'); token=regexprep(token,'_+','_');
    token=regexprep(token,'^[._-]+|[._-]+$',''); if isempty(token), token='run'; end
end

function [outDir,runId]=make_unique_run_directory(parentDir,requestedRunId)
    if ~exist(parentDir,'dir'), mkdir(parentDir); end
    runId=requestedRunId; outDir=fullfile(parentDir,runId); suffix=1;
    while exist(outDir,'dir')
        suffix=suffix+1; runId=sprintf('%s_%02d',requestedRunId,suffix); outDir=fullfile(parentDir,runId);
    end
    [ok,msg]=mkdir(outDir);
    if ~ok, error('diagnose_fast_ami_high_snr:CreateDirFailed','Could not create %s: %s',outDir,msg); end
end

function token=num_token(v,nDec)
    token=sprintf(['%0.' num2str(nDec) 'f'],double(v)); token=strrep(token,'-','m'); token=strrep(token,'.','p');
end

function files=save_figure_pair(fig,figDir,baseName)
    if ~exist(figDir,'dir'), mkdir(figDir); end
    pngPath=fullfile(figDir,[baseName '.png']); figPath=fullfile(figDir,[baseName '.fig']);
    exportgraphics(fig,pngPath,'Resolution',220); savefig(fig,figPath); files={pngPath,figPath};
end
