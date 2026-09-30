function bundle = diagnose_production_dt_policy_accuracy(varargin)
%DIAGNOSE_PRODUCTION_DT_POLICY_ACCURACY  Commit-I.9 full-grid policy sweep.
%
% Repeats the deterministic no-SA production-grid coverage check using the
% Commit-I.9 per-channel-case log-h spacing policy.  The candidate-adaptive
% y-grid remains fixed at scale=1.1.  For every geometry this function also
% evaluates the old all-cases dt=0.01 setting, so the effect of the policy is
% visible directly against the same independent validator.
%
% Default production grid:
%   M={4,8,16,32}, SNR=5:5:30 dB, sigma_X^2={0,0.1,0.2,0.3}
%   P_avg=1, d_min=0.01, pinZero=true
%   adaptive y-grid scale=1.1
%   dt policy: base 0.01, refined 0.005
%
% As in Commit I.7, three deterministic feasible geometries are exercised:
% uniform, upper-quarter, and upper-eighth. Duplicate geometries at small M
% are removed, giving 264 rows on the default grid.
%
% No SA is run. dt is selected once per physical channel point and is reused
% for every geometry at that point.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'MVec',[4 8 16 32],@positive_integer_vector);
    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'SNRVec',5:5:30,@finite_vector);
    addParameter(p,'TurbulenceVec',[0 0.1 0.2 0.3],@nonnegative_vector);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'pinZero',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'LoghDtPolicy','case-adaptive',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'BaseDt',0.01,@positive_scalar);
    addParameter(p,'RefinedDt',0.005,@positive_scalar);
    addParameter(p,'FixedDt',0.01,@positive_scalar);
    addParameter(p,'loghSpanSigma',8,@positive_scalar);
    addParameter(p,'loghYBlockSize',512,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) positive_scalar(v)&&v>=1);
    addParameter(p,'AbsErrorTolerance',1e-3,@positive_scalar);
    addParameter(p,'GeometryModes',{'uniform','upper-quarter','upper-eighth'},@geometry_modes_valid);
    addParameter(p,'UseParallelCases',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'WriteCSV',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:}); o=p.Results;

    if ~o.pinZero
        error('diagnose_production_dt_policy_accuracy:PinZeroRequired', ...
            'Commit-I.9 deterministic stress geometries are defined for pinZero=true.');
    end

    MVec=unique(double(o.MVec(:).'),'stable');
    snrVec=unique(double(o.SNRVec(:).'),'stable');
    sigVec=unique(double(o.TurbulenceVec(:).'),'stable');
    modes=cellstr(string(o.GeometryModes(:).'));
    Pavg=double(o.P_avg); dmin=double(o.minGap);
    validate_spacing(MVec,Pavg,dmin);

    [iM,iSig,iSNR]=ndgrid(1:numel(MVec),1:numel(sigVec),1:numel(snrVec));
    tM=iM(:); tSig=iSig(:); tSNR=iSNR(:); nPhysical=numel(tM);

    % Select dt once per physical channel point before any constellation is scored.
    dtGrid=nan(numel(MVec),numel(sigVec),numel(snrVec));
    dtMeta=cell(size(dtGrid));
    for a=1:numel(MVec)
        for b=1:numel(sigVec)
            for c=1:numel(snrVec)
                [dtGrid(a,b,c),dtMeta{a,b,c}]=select_logh_dt( ...
                    MVec(a),snrVec(c),sigVec(b), ...
                    'Policy',o.LoghDtPolicy,'BaseDt',o.BaseDt, ...
                    'RefinedDt',o.RefinedDt,'FixedDt',o.FixedDt);
            end
        end
    end

    resultsRoot=char(o.ResultsRoot);
    parentDir=fullfile(resultsRoot,'production_dt_policy_accuracy');
    label=sanitize_label(o.RunLabel);
    token=sprintf('policy%s_base%s_ref%s_scale%s_tol%s', ...
        sanitize_label(o.LoghDtPolicy),num_token(o.BaseDt,5),num_token(o.RefinedDt,5), ...
        num_token(o.AdaptiveScale,2),num_token(o.AbsErrorTolerance,5));
    utcToken=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId]=make_unique_dir(parentDir,sprintf('%s_%s_%s',label,token,utcToken));

    runMeta=fso_result_utils.run_metadata(mfilename);
    runMeta.runId=runId; runMeta.runLabel=char(string(o.RunLabel)); runMeta.outputDirectory=outDir;

    refinedPhysical=sum(cellfun(@(m) m.isRefined,dtMeta(:)));
    fprintf('\n============================================================\n');
    fprintf('Commit I.9 - dynamic log-h dt production accuracy sweep\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | Pavg=%.4g | d_min=%.4g\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),Pavg,dmin);
    fprintf('dt policy=%s | base=%.6g | refined=%.6g | refined physical cases=%d/%d\n', ...
        char(string(o.LoghDtPolicy)),o.BaseDt,o.RefinedDt,refinedPhysical,nPhysical);
    fprintf('candidate-adaptive y-grid scale=%.4g | tolerance=%.3e bit/symbol\n', ...
        o.AdaptiveScale,o.AbsErrorTolerance);
    fprintf('Geometry modes=%s | no SA is run\n',strjoin(modes,', '));
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    caseRows=cell(nPhysical,1);
    doPar=logical(o.UseParallelCases);
    if doPar
        try, doPar=license('test','Distrib_Computing_Toolbox')~=0; catch, doPar=false; end
    end
    if doPar
        pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes'); end
        fprintf('Outer physical-point parallelism: %d process workers\n',pool.NumWorkers);
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~) tick());
        parfor k=1:nPhysical
            caseRows{k}=evaluate_point(MVec(tM(k)),snrVec(tSNR(k)),sigVec(tSig(k)), ...
                Pavg,dmin,modes,dtGrid(tM(k),tSig(k),tSNR(k)), ...
                dtMeta{tM(k),tSig(k),tSNR(k)},o);
            send(Q,1);
        end
    else
        fprintf('Outer physical-point parallelism: disabled\n');
        for k=1:nPhysical
            caseRows{k}=evaluate_point(MVec(tM(k)),snrVec(tSNR(k)),sigVec(tSig(k)), ...
                Pavg,dmin,modes,dtGrid(tM(k),tSig(k),tSNR(k)), ...
                dtMeta{tM(k),tSig(k),tSNR(k)},o);
            fprintf('  completed %d/%d\n',k,nPhysical);
        end
    end

    T=vertcat(caseRows{:});
    summary=summarize_table(T,double(o.AbsErrorTolerance),refinedPhysical,nPhysical);

    bundle=struct();
    bundle.meta=runMeta;
    bundle.experiment='production_dt_policy_accuracy';
    bundle.run=struct('id',runId,'label',char(string(o.RunLabel)),'outputDirectory',outDir);
    bundle.axes=struct('M',MVec,'P_avg',Pavg,'SNR_dB',snrVec, ...
        'sigma_X_sq',sigVec,'geometryModes',{modes});
    bundle.policy=struct('name',char(string(o.LoghDtPolicy)),'ruleVersion','I9-v1', ...
        'baseDt',double(o.BaseDt),'refinedDt',double(o.RefinedDt),'fixedDt',double(o.FixedDt), ...
        'dtGrid',dtGrid,'dtMeta',{dtMeta},'refinedPhysicalCases',refinedPhysical);
    bundle.method=struct('fastFadingMethod','logh','loghSpanSigma',double(o.loghSpanSigma), ...
        'loghYBlockSize',double(o.loghYBlockSize),'yGridMode','candidate-adaptive', ...
        'adaptiveScale',double(o.AdaptiveScale),'absErrorTolerance',double(o.AbsErrorTolerance));
    bundle.results=T; bundle.summary=summary;

    matPath=fullfile(outDir,'productionDtPolicyAccuracy.mat');
    save(matPath,'bundle','-v7.3');
    if o.WriteCSV, writetable(T,fullfile(outDir,'productionDtPolicyAccuracy.csv')); end
    fprintf('\nSaved I.9 policy accuracy bundle: %s\n',matPath);
    print_summary(summary,T);

    function tick()
        doneCount=doneCount+1; fprintf('  completed %d/%d\n',doneCount,nPhysical);
    end
end


function T=evaluate_point(M,SNRdB,sigmaX2,Pavg,dmin,modes,selectedDt,dtMeta,o)
    xMaxBound=fso_result_utils.feasible_xmax_bound(M,Pavg,dmin,true);
    cfg=build_fso_config(M,Pavg,SNRdB,sigmaX2, ...
        'FastFadingMethod','logh','loghDt',selectedDt, ...
        'loghSpanSigma',o.loghSpanSigma,'loghYBlockSize',o.loghYBlockSize, ...
        'YGridMode','candidate-adaptive','YGridScale',o.AdaptiveScale, ...
        'xMaxBound',xMaxBound,'saMaxIter',1,'saNStarts',1,'saUseParallel',false, ...
        'minGap',dmin,'pinZero',true,'seedInit',1,'logEvery',0,'historyEvery',1);

    [X,names]=deterministic_geometries(M,Pavg,dmin,modes); nG=numel(names);
    Mcol=repmat(M,nG,1); SigmaX2=repmat(sigmaX2,nG,1); SNRdBcol=repmat(SNRdB,nG,1);
    Geometry=strings(nG,1); SelectedDt=repmat(selectedDt,nG,1); IsRefined=repmat(dtMeta.isRefined,nG,1);
    DtReason=repmat(string(dtMeta.reason),nG,1);
    XMax=nan(nG,1); XMaxToBound=nan(nG,1); ValidatorAMI=nan(nG,1);
    BaseAMI=nan(nG,1); PolicyAMI=nan(nG,1); BaseError=nan(nG,1); PolicyError=nan(nG,1);
    PolicyMinusBase=nan(nG,1); ValidatorSeconds=nan(nG,1); BaseSeconds=nan(nG,1); PolicySeconds=nan(nG,1);
    YGridN=nan(nG,1); AdaptiveSupport=nan(nG,1);

    for g=1:nG
        x=X{g}; Geometry(g)=string(names{g});
        AMI_functions.assert_constellation_feasible(x,cfg,1e-10);
        XMax(g)=max(abs(x)); XMaxToBound(g)=XMax(g)/xMaxBound;
        tv=tic; ValidatorAMI(g)=cfg.AMI_Validator(x(:)); ValidatorSeconds(g)=toc(tv);

        % Build the same candidate-adaptive y-grid once, then isolate dt only.
        yGrid=AMI_functions.build_noCSI_y_grid(cfg,o.AdaptiveScale*XMax(g));
        YGridN(g)=numel(yGrid); AdaptiveSupport(g)=o.AdaptiveScale*XMax(g);

        tb=tic;
        BaseAMI(g)=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,o.BaseDt, ...
            o.loghSpanSigma,yGrid,o.loghYBlockSize);
        BaseSeconds(g)=toc(tb);

        tp=tic;
        PolicyAMI(g)=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,selectedDt, ...
            o.loghSpanSigma,yGrid,o.loghYBlockSize);
        PolicySeconds(g)=toc(tp);

        BaseError(g)=BaseAMI(g)-ValidatorAMI(g);
        PolicyError(g)=PolicyAMI(g)-ValidatorAMI(g);
        PolicyMinusBase(g)=PolicyAMI(g)-BaseAMI(g);
    end

    T=table(Mcol,SigmaX2,SNRdBcol,Geometry,SelectedDt,IsRefined,DtReason, ...
        XMax,XMaxToBound,ValidatorAMI,BaseAMI,PolicyAMI,BaseError,PolicyError, ...
        PolicyMinusBase,ValidatorSeconds,BaseSeconds,PolicySeconds,YGridN,AdaptiveSupport, ...
        'VariableNames',{'M','SigmaX2','SNRdB','Geometry','SelectedDt','IsRefined','DtReason', ...
        'XMax','XMaxToBound','ValidatorAMI','BaseAMI','PolicyAMI','BaseError','PolicyError', ...
        'PolicyMinusBase','ValidatorSeconds','BaseSeconds','PolicySeconds','YGridN','AdaptiveSupport'});
end


function [X,names]=deterministic_geometries(M,Pavg,dmin,modes)
    X={}; names={};
    for j=1:numel(modes)
        mode=lower(string(modes{j}));
        switch mode
            case "uniform", x=define_constellation(M,Pavg,true,"mean");
            case "upper-quarter", x=clustered_geometry(M,Pavg,dmin,max(1,ceil(M/4)));
            case "upper-eighth", x=clustered_geometry(M,Pavg,dmin,max(1,ceil(M/8)));
            case "upper-one", x=clustered_geometry(M,Pavg,dmin,1);
            otherwise, error('diagnose_production_dt_policy_accuracy:BadGeometryMode','Unknown geometry mode: %s',mode);
        end
        x=real(x(:)); duplicate=false;
        for k=1:numel(X)
            if numel(X{k})==numel(x) && max(abs(X{k}-x))<=1e-12*max(1,max(abs(x)))
                duplicate=true; break;
            end
        end
        if ~duplicate, X{end+1,1}=x; names{end+1,1}=char(mode); end %#ok<AGROW>
    end
end

function x=clustered_geometry(M,Pavg,dmin,K)
    K=min(max(1,round(K)),M-1); b=(0:M-1).'*dmin; surplus=Pavg-mean(b);
    if surplus < -100*eps(max(1,Pavg)), error('diagnose_production_dt_policy_accuracy:InfeasibleGap','Minimum-gap skeleton exceeds Pavg.'); end
    surplus=max(0,surplus); q=zeros(M,1); q(end-K+1:end)=M*surplus/K; x=b+q; x(1)=0;
    err=Pavg-mean(x); if abs(err)>0, x(end-K+1:end)=x(end-K+1:end)+M*err/K; end
end

function S=summarize_table(T,tol,nRefinedPhysical,nPhysical)
    absB=abs(T.BaseError); absP=abs(T.PolicyError);
    [S.maxAbsBaseError,idxB]=max(absB,[],'omitnan');
    [S.maxAbsPolicyError,idxP]=max(absP,[],'omitnan');
    S.meanAbsBaseError=mean(absB,'omitnan'); S.meanAbsPolicyError=mean(absP,'omitnan');
    S.maxAbsPolicyMinusBase=max(abs(T.PolicyMinusBase),[],'omitnan');
    S.absErrorTolerance=tol; S.baseAllPass=all(absB<=tol|isnan(absB)); S.policyAllPass=all(absP<=tol|isnan(absP));
    S.nRows=height(T); S.nBaseFailures=sum(absB>tol); S.nPolicyFailures=sum(absP>tol);
    S.nRefinedRows=sum(T.IsRefined); S.nRefinedPhysicalCases=nRefinedPhysical; S.nPhysicalCases=nPhysical;
    S.meanBaseSeconds=mean(T.BaseSeconds,'omitnan'); S.meanPolicySeconds=mean(T.PolicySeconds,'omitnan');
    S.totalBaseSeconds=sum(T.BaseSeconds,'omitnan'); S.totalPolicySeconds=sum(T.PolicySeconds,'omitnan');
    S.policyRuntimeRatioToBase=S.totalPolicySeconds/max(S.totalBaseSeconds,eps);
    S.worstBaseRow=idxB; S.worstPolicyRow=idxP;
end

function print_summary(S,T)
    fprintf('\nI.9 DYNAMIC-DT PRODUCTION ACCURACY SUMMARY\n');
    fprintf('Rows evaluated                    : %d\n',S.nRows);
    fprintf('Physical cases refined            : %d/%d\n',S.nRefinedPhysicalCases,S.nPhysicalCases);
    fprintf('Tolerance                         : %.3e bit/symbol\n',S.absErrorTolerance);
    fprintf('All-dt=0.01 mean/max |error|      : %.3e / %.3e | failures=%d\n',S.meanAbsBaseError,S.maxAbsBaseError,S.nBaseFailures);
    fprintf('Dynamic policy mean/max |error|   : %.3e / %.3e | failures=%d\n',S.meanAbsPolicyError,S.maxAbsPolicyError,S.nPolicyFailures);
    fprintf('Dynamic policy all rows pass      : %d\n',S.policyAllPass);
    fprintf('Max |policy-base AMI|             : %.3e\n',S.maxAbsPolicyMinusBase);
    fprintf('Total policy/base evaluator time  : %.3fx\n',S.policyRuntimeRatioToBase);
    if ~isempty(S.worstPolicyRow) && isfinite(S.worstPolicyRow)
        r=T(S.worstPolicyRow,:);
        fprintf(['Worst policy case: M=%d SNR=%.1f sigma_X^2=%.3f geometry=%s ' ...
            '| dt=%.5g | error=%+.3e | refined=%d\n'], ...
            r.M,r.SNRdB,r.SigmaX2,char(r.Geometry),r.SelectedDt,r.PolicyError,r.IsRefined);
    end
end

function validate_spacing(MVec,Pavg,dmin)
    req=(MVec-1)*dmin/2; bad=find(req>Pavg+100*eps(max(1,Pavg)),1,'first');
    if ~isempty(bad), error('diagnose_production_dt_policy_accuracy:InfeasibleGap', ...
        'M=%d and d_min=%.6g require mean power %.6g > P_avg=%.6g.',MVec(bad),dmin,req(bad),Pavg); end
end
function tf=geometry_modes_valid(v)
    if isstring(v), v=cellstr(v); end; tf=iscell(v)&&~isempty(v); if ~tf, return; end
    allowed=["uniform","upper-quarter","upper-eighth","upper-one"];
    try, tf=all(ismember(lower(string(v)),allowed)); catch, tf=false; end
end
function t=sanitize_label(v), t=strtrim(char(string(v))); if isempty(t),t='run';return;end; t=regexprep(t,'[^A-Za-z0-9._-]+','_'); t=regexprep(t,'_+','_'); t=regexprep(t,'^[._-]+|[._-]+$',''); if isempty(t),t='run';end; end
function [d,id]=make_unique_dir(parent,id), if ~exist(parent,'dir'),mkdir(parent);end; d=fullfile(parent,id);base=id;n=1;while exist(d,'dir'),n=n+1;id=sprintf('%s_%02d',base,n);d=fullfile(parent,id);end;[ok,msg]=mkdir(d);if ~ok,error('diagnose_production_dt_policy_accuracy:mkdir','%s',msg);end;end
function token=num_token(v,nDec), token=sprintf(['%0.' num2str(nDec) 'f'],double(v));token=strrep(token,'-','m');token=strrep(token,'+','p');token=strrep(token,'.','p');end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=nonnegative_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;end
function tf=positive_integer(v),tf=positive_scalar(v)&&mod(v,1)==0;end
function tf=positive_integer_vector(v),tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=1)&&all(mod(v,1)==0);end
function tf=finite_vector(v),tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v));end
function tf=nonnegative_vector(v),tf=finite_vector(v)&&all(v>=0);end
