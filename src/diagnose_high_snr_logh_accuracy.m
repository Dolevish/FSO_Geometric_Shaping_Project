function bundle = diagnose_high_snr_logh_accuracy(varargin)
%DIAGNOSE_HIGH_SNR_LOGH_ACCURACY  Commit-I.8 high-SNR integration diagnostic.
%
% Isolates the residual fast-AMI error exposed by the Commit-I.7 full-grid
% sweep. No SA is run and no production default is changed.
%
% The diagnostic independently refines:
%   (1) the fading integration step dt in t=ln(h), and
%   (2) the existing nonuniform outer y-grid by an integer refinement factor.
%
% Both fixed-bound and candidate-adaptive y-support modes are tested against
% exactly the same independent validator. This distinguishes a log-h
% discretization error from an outer-y discretization/support error.
%
% Default stress cases are the I.7 high-SNR cases at/near the 1e-3 target:
%   [M, SNR_dB, sigma_X^2] =
%       [32, 30, 0.30]
%       [32, 30, 0.20]
%       [16, 30, 0.30]
%       [32, 30, 0.10]
%
% Name-value options:
%   'Cases'             N-by-3 numeric matrix [M,SNR_dB,sigma_X^2]
%   'GeometryMode'      'upper-eighth' (default), 'upper-quarter',
%                       'upper-one', or 'uniform'
%   'P_avg'             average optical intensity (default 1)
%   'minGap'            minimum adjacent level spacing (default 0.01)
%   'DtVec'             log-h steps (default [0.01 .005 .0025 .00125])
%   'YRefineVec'        integer subdivision factors (default [1 2 4])
%   'TSpanSigma'        log-h half-span in sigma_t (default 8)
%   'YBlockSize'        y processing block size (default 512)
%   'AdaptiveScale'     adaptive support multiplier (default 1.1)
%   'AbsErrorTolerance' absolute error target (default 1e-3 bit/symbol)
%   'UseParallelCases'  parallelize independent stress cases (default true)
%   'WriteCSV'          save flat result/summary CSV files (default true)
%   'RunLabel'          optional run label
%   'ResultsRoot'       output root
%
% The y-refinement operation preserves the original y-grid endpoints and
% inserts equally spaced points inside every existing interval. Therefore a
% y-refinement comparison changes resolution only, not support.

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p,'Cases',[32 30 .30; 32 30 .20; 16 30 .30; 32 30 .10],@cases_valid);
    addParameter(p,'GeometryMode','upper-eighth',@geometry_mode_valid);
    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'DtVec',[0.01 0.005 0.0025 0.00125],@(v) numeric_vector(v)&&all(v>0));
    addParameter(p,'YRefineVec',[1 2 4],@(v) numeric_vector(v)&&all(v>=1)&&all(mod(v,1)==0));
    addParameter(p,'TSpanSigma',8,@positive_scalar);
    addParameter(p,'YBlockSize',512,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) positive_scalar(v)&&v>=1);
    addParameter(p,'AbsErrorTolerance',1e-3,@positive_scalar);
    addParameter(p,'UseParallelCases',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'WriteCSV',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:});
    o = p.Results;

    cases = double(o.Cases);
    dtVec = unique(double(o.DtVec(:).'),'stable');
    yRefineVec = unique(double(o.YRefineVec(:).'),'stable');
    geometryMode = char(lower(string(o.GeometryMode)));
    Pavg = double(o.P_avg);
    dmin = double(o.minGap);
    tol = double(o.AbsErrorTolerance);

    for k=1:size(cases,1)
        validate_case(cases(k,:),Pavg,dmin);
    end

    parentDir = fullfile(char(o.ResultsRoot),'high_snr_logh_accuracy');
    label = sanitize_label(o.RunLabel);
    runToken = sprintf('DT%s_Yref%s_scale%s_tol%s', ...
        join_decimals(dtVec),join_ints(yRefineVec),num_token(o.AdaptiveScale,2),num_token(tol,5));
    utcToken = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId] = make_unique_dir(parentDir,sprintf('%s_%s_%s',label,runToken,utcToken));

    runMeta = fso_result_utils.run_metadata(mfilename);
    runMeta.runId = runId;
    runMeta.runLabel = char(string(o.RunLabel));
    runMeta.outputDirectory = outDir;

    fprintf('\n============================================================\n');
    fprintf('Commit I.8 - high-SNR log-h / y-grid accuracy diagnostic\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('Cases [M SNR sigma_X^2]:\n'); disp(cases);
    fprintf('Geometry: %s | Pavg=%.4g | d_min=%.4g\n',geometryMode,Pavg,dmin);
    fprintf('log-h dt values: %s | span=+/-%.3g sigma_t\n',mat2str(dtVec),o.TSpanSigma);
    fprintf('y-grid refinement factors: %s | adaptive scale=%.4g\n',mat2str(yRefineVec),o.AdaptiveScale);
    fprintf('Tolerance: %.3e bit/symbol | no SA is run\n',tol);
    fprintf('Each y-refinement preserves its base support endpoints.\n');
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    nCases = size(cases,1);
    caseTables = cell(nCases,1);
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
        afterEach(Q,@(~) tick());
        parfor k=1:nCases
            caseTables{k} = evaluate_case(cases(k,:),geometryMode,Pavg,dmin,dtVec,yRefineVec,o);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k=1:nCases
            caseTables{k} = evaluate_case(cases(k,:),geometryMode,Pavg,dmin,dtVec,yRefineVec,o);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    T = vertcat(caseTables{:});
    S = build_summary(T,dtVec,yRefineVec,tol);

    bundle = struct();
    bundle.meta = runMeta;
    bundle.experiment = 'high_snr_logh_accuracy';
    bundle.run = struct('id',runId,'label',char(string(o.RunLabel)),'outputDirectory',outDir);
    bundle.cases = cases;
    bundle.settings = struct('geometryMode',geometryMode,'P_avg',Pavg,'minGap',dmin, ...
        'dtVec',dtVec,'yRefineVec',yRefineVec,'TSpanSigma',double(o.TSpanSigma), ...
        'YBlockSize',double(o.YBlockSize),'adaptiveScale',double(o.AdaptiveScale), ...
        'absErrorTolerance',tol,'outerParallelism',doPar);
    bundle.results = T;
    bundle.summary = S;

    matPath = fullfile(outDir,'highSNRLoghAccuracy.mat');
    save(matPath,'bundle','-v7.3');
    if o.WriteCSV
        writetable(T,fullfile(outDir,'highSNRLoghAccuracy.csv'));
        writetable(S.bySetting,fullfile(outDir,'highSNRLoghAccuracy_summary.csv'));
    end

    fprintf('\nSaved I.8 diagnostic bundle: %s\n',matPath);
    print_summary(S,T);

    function tick()
        doneCount = doneCount + 1;
        fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function T = evaluate_case(caseRow,geometryMode,Pavg,dmin,dtVec,yRefineVec,o)
    M = caseRow(1); SNRdB = caseRow(2); sigmaX2 = caseRow(3);
    x = make_geometry(M,Pavg,dmin,geometryMode);
    xMaxBound = fso_result_utils.feasible_xmax_bound(M,Pavg,dmin,true);

    cfg = build_fso_config(M,Pavg,SNRdB,sigmaX2, ...
        'FastFadingMethod','logh','loghDt',dtVec(1), ...
        'loghSpanSigma',o.TSpanSigma,'loghYBlockSize',o.YBlockSize, ...
        'YGridMode','fixed-bound','YGridScale',o.AdaptiveScale, ...
        'xMaxBound',xMaxBound,'saMaxIter',1,'saNStarts',1,'saUseParallel',false, ...
        'minGap',dmin,'pinZero',true,'seedInit',1,'logEvery',0,'historyEvery',1);
    AMI_functions.assert_constellation_feasible(x,cfg,1e-10);

    tv=tic;
    validatorAMI = cfg.AMI_Validator(x(:));
    validatorSeconds = toc(tv);

    fixedBase = cfg.y_grid(:);
    adaptiveSupport = double(o.AdaptiveScale)*max(abs(x));
    adaptiveBase = AMI_functions.build_noCSI_y_grid(cfg,adaptiveSupport);

    nRows = 2*numel(dtVec)*numel(yRefineVec);
    Mcol = repmat(M,nRows,1);
    SigmaX2 = repmat(sigmaX2,nRows,1);
    SNRcol = repmat(SNRdB,nRows,1);
    Geometry = repmat(string(geometryMode),nRows,1);
    Mode = strings(nRows,1);
    Dt = nan(nRows,1);
    YRefine = nan(nRows,1);
    XMax = repmat(max(abs(x)),nRows,1);
    XMaxBound = repmat(xMaxBound,nRows,1);
    XMaxToBound = repmat(max(abs(x))/xMaxBound,nRows,1);
    YSupportMax = nan(nRows,1);
    YGridN = nan(nRows,1);
    ValidatorAMI = repmat(validatorAMI,nRows,1);
    FastAMI = nan(nRows,1);
    SignedError = nan(nRows,1);
    AbsError = nan(nRows,1);
    EvalSeconds = nan(nRows,1);
    ValidatorSeconds = repmat(validatorSeconds,nRows,1);

    row=0;
    for im=1:2
        if im==1
            mode="fixed-bound"; baseGrid=fixedBase;
        else
            mode="candidate-adaptive"; baseGrid=adaptiveBase;
        end
        for iy=1:numel(yRefineVec)
            yr = refine_grid(baseGrid,yRefineVec(iy));
            for id=1:numel(dtVec)
                row=row+1;
                tt=tic;
                v=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,dtVec(id),o.TSpanSigma,yr,o.YBlockSize);
                EvalSeconds(row)=toc(tt);
                Mode(row)=mode;
                Dt(row)=dtVec(id);
                YRefine(row)=yRefineVec(iy);
                YSupportMax(row)=yr(end);
                YGridN(row)=numel(yr);
                FastAMI(row)=v;
                SignedError(row)=v-validatorAMI;
                AbsError(row)=abs(SignedError(row));
            end
        end
    end

    T=table(Mcol,SigmaX2,SNRcol,Geometry,Mode,Dt,YRefine,XMax,XMaxBound,XMaxToBound, ...
        YSupportMax,YGridN,ValidatorAMI,FastAMI,SignedError,AbsError,EvalSeconds,ValidatorSeconds, ...
        'VariableNames',{'M','SigmaX2','SNRdB','Geometry','Mode','Dt','YRefine','XMax', ...
        'XMaxBound','XMaxToBound','YSupportMax','YGridN','ValidatorAMI','FastAMI', ...
        'SignedError','AbsError','EvalSeconds','ValidatorSeconds'});
end


function S = build_summary(T,dtVec,yRefineVec,tol)
    modes=["fixed-bound","candidate-adaptive"];
    n=2*numel(dtVec)*numel(yRefineVec);
    Mode=strings(n,1); Dt=nan(n,1); YRefine=nan(n,1);
    MeanAbsError=nan(n,1); MaxAbsError=nan(n,1); MedianSeconds=nan(n,1);
    MeanGridN=nan(n,1); AllPass=false(n,1); nFail=nan(n,1);
    row=0;
    for im=1:2
        for iy=1:numel(yRefineVec)
            for id=1:numel(dtVec)
                row=row+1;
                mask=T.Mode==modes(im) & abs(T.Dt-dtVec(id))<1e-15 & T.YRefine==yRefineVec(iy);
                a=T.AbsError(mask);
                Mode(row)=modes(im); Dt(row)=dtVec(id); YRefine(row)=yRefineVec(iy);
                MeanAbsError(row)=mean(a,'omitnan'); MaxAbsError(row)=max(a,[],'omitnan');
                MedianSeconds(row)=median(T.EvalSeconds(mask),'omitnan');
                MeanGridN(row)=mean(T.YGridN(mask),'omitnan');
                AllPass(row)=all(a<=tol); nFail(row)=sum(a>tol);
            end
        end
    end
    S.bySetting=table(Mode,Dt,YRefine,MeanAbsError,MaxAbsError,MedianSeconds,MeanGridN,AllPass,nFail);
    S.absErrorTolerance=tol;

    baseMask=T.Mode=="candidate-adaptive" & T.YRefine==1;
    passDt=[];
    for d=dtVec
        a=T.AbsError(baseMask & abs(T.Dt-d)<1e-15);
        if ~isempty(a) && all(a<=tol), passDt(end+1)=d; end %#ok<AGROW>
    end
    if isempty(passDt), S.coarsestPassingAdaptiveDtAtBaseY=NaN;
    else, S.coarsestPassingAdaptiveDtAtBaseY=max(passDt); end

    d0=max(dtVec);
    passY=[];
    for r=yRefineVec
        a=T.AbsError(T.Mode=="candidate-adaptive" & abs(T.Dt-d0)<1e-15 & T.YRefine==r);
        if ~isempty(a) && all(a<=tol), passY(end+1)=r; end %#ok<AGROW>
    end
    if isempty(passY), S.smallestPassingAdaptiveYRefineAtCoarsestDt=NaN;
    else, S.smallestPassingAdaptiveYRefineAtCoarsestDt=min(passY); end

    finestDt=min(dtVec); maxRef=max(yRefineVec);
    S.maxErrorAdaptiveCurrent = max(T.AbsError(T.Mode=="candidate-adaptive" & abs(T.Dt-d0)<1e-15 & T.YRefine==1),[],'omitnan');
    S.maxErrorAdaptiveFinestDtBaseY = max(T.AbsError(T.Mode=="candidate-adaptive" & abs(T.Dt-finestDt)<1e-15 & T.YRefine==1),[],'omitnan');
    S.maxErrorAdaptiveCoarsestDtMaxY = max(T.AbsError(T.Mode=="candidate-adaptive" & abs(T.Dt-d0)<1e-15 & T.YRefine==maxRef),[],'omitnan');
    S.maxErrorAdaptiveFinestBoth = max(T.AbsError(T.Mode=="candidate-adaptive" & abs(T.Dt-finestDt)<1e-15 & T.YRefine==maxRef),[],'omitnan');
    S.nCases=numel(unique(string(T.M)+"_"+string(T.SNRdB)+"_"+string(T.SigmaX2)));
end


function print_summary(S,T)
    fprintf('\nI.8 HIGH-SNR ACCURACY SUMMARY\n');
    fprintf('Stress cases: %d | tolerance: %.3e bit/symbol\n',S.nCases,S.absErrorTolerance);
    fprintf('Adaptive current (coarsest dt, yRef=1) max error : %.3e\n',S.maxErrorAdaptiveCurrent);
    fprintf('Adaptive finest dt, yRef=1 max error             : %.3e\n',S.maxErrorAdaptiveFinestDtBaseY);
    fprintf('Adaptive coarsest dt, max y-refine max error      : %.3e\n',S.maxErrorAdaptiveCoarsestDtMaxY);
    fprintf('Adaptive finest dt + max y-refine max error       : %.3e\n',S.maxErrorAdaptiveFinestBoth);
    fprintf('Coarsest adaptive dt passing all cases at yRef=1  : %.6g\n',S.coarsestPassingAdaptiveDtAtBaseY);
    fprintf('Smallest y-refine passing at coarsest dt          : %.6g\n',S.smallestPassingAdaptiveYRefineAtCoarsestDt);

    fprintf('\nAdaptive yRef=1 by dt:\n');
    A=S.bySetting(S.bySetting.Mode=="candidate-adaptive" & S.bySetting.YRefine==1,:);
    A=sortrows(A,'Dt','descend');
    disp(A(:,{'Dt','MeanAbsError','MaxAbsError','MedianSeconds','AllPass','nFail'}));

    fprintf('Adaptive coarsest dt by y refinement:\n');
    d0=max(S.bySetting.Dt);
    B=S.bySetting(S.bySetting.Mode=="candidate-adaptive" & abs(S.bySetting.Dt-d0)<1e-15,:);
    B=sortrows(B,'YRefine','ascend');
    disp(B(:,{'YRefine','MeanAbsError','MaxAbsError','MedianSeconds','MeanGridN','AllPass','nFail'}));

    [~,idx]=sort(T.AbsError,'descend');
    nShow=min(8,numel(idx));
    fprintf('Worst %d individual evaluations:\n',nShow);
    disp(T(idx(1:nShow),{'M','SigmaX2','SNRdB','Mode','Dt','YRefine','AbsError','YGridN','EvalSeconds'}));
end


function y2=refine_grid(y,factor)
    y=double(y(:)); factor=double(factor);
    if factor==1, y2=y; return; end
    dy=diff(y);
    alpha=(0:factor-1)/factor;
    a=y(1:end-1)+dy.*alpha;
    y2=[reshape(a.',[],1); y(end)];
    if any(diff(y2)<=0)
        error('diagnose_high_snr_logh_accuracy:RefineGrid','Refined y-grid is not strictly increasing.');
    end
end


function x=make_geometry(M,Pavg,dmin,mode)
    mode=lower(string(mode));
    switch mode
        case "uniform"
            x=define_constellation(M,Pavg,true,"mean");
        case "upper-quarter"
            x=clustered_geometry(M,Pavg,dmin,max(1,ceil(M/4)));
        case "upper-eighth"
            x=clustered_geometry(M,Pavg,dmin,max(1,ceil(M/8)));
        case "upper-one"
            x=clustered_geometry(M,Pavg,dmin,1);
        otherwise
            error('diagnose_high_snr_logh_accuracy:BadGeometry','Unknown geometry mode %s.',mode);
    end
    x=real(x(:));
end

function x=clustered_geometry(M,Pavg,dmin,K)
    K=min(max(1,round(K)),M-1);
    b=(0:M-1).'*dmin;
    surplus=Pavg-mean(b);
    if surplus < -100*eps(max(1,Pavg))
        error('diagnose_high_snr_logh_accuracy:InfeasibleGap','Minimum-gap skeleton exceeds P_avg.');
    end
    q=zeros(M,1); q(end-K+1:end)=M*max(0,surplus)/K;
    x=b+q; x(1)=0;
    err=Pavg-mean(x);
    x(end-K+1:end)=x(end-K+1:end)+M*err/K;
end

function validate_case(c,Pavg,dmin)
    M=c(1); snr=c(2); sig=c(3);
    if ~(isfinite(M)&&M>=2&&mod(M,1)==0&&isfinite(snr)&&isfinite(sig)&&sig>=0)
        error('diagnose_high_snr_logh_accuracy:BadCase','Each case must be [integer M>=2, finite SNR, sigma_X^2>=0].');
    end
    req=(M-1)*dmin/2;
    if req>Pavg+100*eps(max(1,Pavg))
        error('diagnose_high_snr_logh_accuracy:InfeasibleGap','M=%d requires minimum mean %.6g > P_avg %.6g.',M,req,Pavg);
    end
end

function tf=cases_valid(v), tf=isnumeric(v)&&ismatrix(v)&&size(v,2)==3&&~isempty(v)&&all(isfinite(v(:))); end
function tf=geometry_mode_valid(v), tf=(ischar(v)||(isstring(v)&&isscalar(v)))&&ismember(lower(string(v)),["uniform","upper-quarter","upper-eighth","upper-one"]); end
function tf=numeric_vector(v), tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v)); end
function tf=positive_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0; end
function tf=nonnegative_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0; end
function tf=positive_integer(v), tf=positive_scalar(v)&&mod(v,1)==0; end

function t=sanitize_label(v)
    t=strtrim(char(string(v))); if isempty(t), t='run'; return; end
    t=regexprep(t,'[^A-Za-z0-9._-]+','_'); t=regexprep(t,'_+','_');
    t=regexprep(t,'^[._-]+|[._-]+$',''); if isempty(t), t='run'; end
end
function [d,id]=make_unique_dir(parent,id)
    if ~exist(parent,'dir'), mkdir(parent); end
    d=fullfile(parent,id); base=id; n=1;
    while exist(d,'dir'), n=n+1; id=sprintf('%s_%02d',base,n); d=fullfile(parent,id); end
    [ok,msg]=mkdir(d); if ~ok, error('diagnose_high_snr_logh_accuracy:mkdir','%s',msg); end
end
function token=num_token(v,nDec)
    token=sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token=strrep(token,'-','m'); token=strrep(token,'+','p'); token=strrep(token,'.','p');
end
function s=join_decimals(v)
    c=cell(1,numel(v)); for k=1:numel(v), c{k}=num_token(v(k),5); end; s=strjoin(c,'-');
end
function s=join_ints(v), s=strjoin(arrayfun(@(x)sprintf('%d',x),v,'UniformOutput',false),'-'); end
