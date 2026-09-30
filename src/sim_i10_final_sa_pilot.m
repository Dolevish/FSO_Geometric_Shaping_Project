function bundle = sim_i10_final_sa_pilot(varargin)
%SIM_I10_FINAL_SA_PILOT  Final pre-production SA reliability pilot.
%
% Commit I.10 is diagnostic-only. It does not modify the I.9 production
% policy. The pilot exercises the five critical high-SNR cases with the
% production driver itself and the Commit-G SA profile/budget:
%   [M SNR_dB sigma_X^2] =
%      [32 30 0.1]
%      [32 30 0.2]
%      [32 30 0.3]
%      [16 30 0.3]
%      [16 30 0.2]   % boundary case that remains at dt=0.01
%
% All five cases run with the I.9 case-adaptive dt policy. In addition, the
% boundary case is rerun with fixed dt=0.005 using the SAME replicate seeds.
% This isolates whether keeping that boundary point at dt=0.01 changes the
% validated SA result in a practically meaningful way.
%
% Default SA budget/profile:
%   4000 iterations x 4 restarts x 2 replicates
%   T0=0.4, Tf=1e-3, BaseStd0=0.15, ItersPerTemp=50
%
% Acceptance checks:
%   - every policy run max restart fast-vs-validator gap <= 1e-3 bit/symbol
%   - no policy run changes the selected restart after validation
%   - boundary policy/reference validated AMI difference <= 1e-3 bit/symbol
%   - paired boundary seeds match exactly
%
% k99 and restart-to-restart validated std are reported, not thresholded.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'Cases',[32 30 0.1;32 30 0.2;32 30 0.3;16 30 0.3;16 30 0.2],@cases_valid);
    addParameter(p,'BoundaryCase',[16 30 0.2],@(v) isnumeric(v)&&isequal(size(v),[1 3])&&all(isfinite(v)));
    addParameter(p,'RunBoundaryReference',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'ReferenceBoundaryDt',0.005,@positive_scalar);

    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'pinZero',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'BaseDt',0.01,@positive_scalar);
    addParameter(p,'RefinedDt',0.005,@positive_scalar);
    addParameter(p,'loghSpanSigma',8,@positive_scalar);
    addParameter(p,'loghYBlockSize',512,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) positive_scalar(v)&&v>=1);

    addParameter(p,'saMaxIter',4000,@positive_integer);
    addParameter(p,'saNStarts',4,@(v) positive_integer(v)&&v>=2);
    addParameter(p,'nReplicates',2,@positive_integer);
    addParameter(p,'T0',0.4,@positive_scalar);
    addParameter(p,'Tf',1e-3,@positive_scalar);
    addParameter(p,'BaseStd0',0.15,@positive_scalar);
    addParameter(p,'ItersPerTemp',50,@positive_integer);
    addParameter(p,'baseSeed',20260928,@nonnegative_scalar);
    addParameter(p,'historyEvery',50,@positive_integer);

    addParameter(p,'FastValidationTolerance',1e-3,@positive_scalar);
    addParameter(p,'BoundaryValidatedDeltaTolerance',1e-3,@positive_scalar);
    addParameter(p,'UseParallelChildRuns',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'WriteCSV',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:}); o=p.Results;

    if ~(o.T0>o.Tf), error('sim_i10_final_sa_pilot:TemperatureOrder','T0 must exceed Tf.'); end
    if o.RefinedDt>o.BaseDt
        error('sim_i10_final_sa_pilot:DtOrder','RefinedDt must not exceed BaseDt.');
    end

    cases=double(o.Cases);
    if size(unique(cases,'rows'),1)~=size(cases,1)
        error('sim_i10_final_sa_pilot:DuplicateCases','Cases must be unique.');
    end
    boundary=double(o.BoundaryCase);
    if o.RunBoundaryReference && ~any(all(abs(cases-boundary)<=1e-12,2))
        error('sim_i10_final_sa_pilot:BoundaryMissing','BoundaryCase must appear in Cases.');
    end

    % Resolve and display the I.9 policy before any SA starts.
    nPolicy=size(cases,1);
    selectedDt=nan(nPolicy,1); selectedMeta=cell(nPolicy,1);
    for k=1:nPolicy
        [selectedDt(k),selectedMeta{k}]=select_logh_dt(cases(k,1),cases(k,2),cases(k,3), ...
            'Policy','case-adaptive','BaseDt',o.BaseDt,'RefinedDt',o.RefinedDt);
    end

    parentDir=fullfile(char(o.ResultsRoot),'i10_final_sa_pilot');
    label=sanitize_label(o.RunLabel);
    token=sprintf('I%d_S%d_R%d_base%s_ref%s',o.saMaxIter,o.saNStarts,o.nReplicates, ...
        num_token(o.BaseDt,4),num_token(o.RefinedDt,4));
    utcToken=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId]=make_unique_dir(parentDir,sprintf('%s_%s_%s',label,token,utcToken));
    childRoot=fullfile(outDir,'production_children'); mkdir(childRoot);

    runMeta=fso_result_utils.run_metadata(mfilename);
    runMeta.runId=runId; runMeta.runLabel=char(string(o.RunLabel)); runMeta.outputDirectory=outDir;

    nTasks=nPolicy+double(o.RunBoundaryReference);
    taskMode=strings(nTasks,1); taskCase=nan(nTasks,3);
    for k=1:nPolicy
        taskMode(k)="policy"; taskCase(k,:)=cases(k,:);
    end
    if o.RunBoundaryReference
        taskMode(end)="boundary-reference"; taskCase(end,:)=boundary;
    end

    fprintf('\n============================================================\n');
    fprintf('Commit I.10 - final pre-production SA pilot\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('Policy cases [M SNR sigma_X^2] and selected dt:\n');
    for k=1:nPolicy
        fprintf('  [%d %.1f %.3f] -> dt=%.6g (%s)\n',cases(k,1),cases(k,2),cases(k,3), ...
            selectedDt(k),ternary(selectedMeta{k}.isRefined,'refined','base'));
    end
    if o.RunBoundaryReference
        fprintf('Boundary reference: [%d %.1f %.3f] with fixed dt=%.6g and identical seeds\n', ...
            boundary(1),boundary(2),boundary(3),o.ReferenceBoundaryDt);
    end
    fprintf('y-grid=candidate-adaptive scale=%.3g | log-h span=+/-%.3g sigma_t\n', ...
        o.AdaptiveScale,o.loghSpanSigma);
    fprintf('SA=%d restarts x %d iterations x %d replicates | T0=%.3g Tf=%.3g step0=%.3g block=%d\n', ...
        o.saNStarts,o.saMaxIter,o.nReplicates,o.T0,o.Tf,o.BaseStd0,o.ItersPerTemp);
    fprintf('Child production runs=%d | total SA optimization cases=%d\n',nTasks,nTasks*o.nReplicates);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    childBundles=cell(nTasks,1);
    doPar=logical(o.UseParallelChildRuns);
    if doPar
        try, doPar=license('test','Distrib_Computing_Toolbox')~=0; catch, doPar=false; end
    end

    if doPar
        pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes'); end
        fprintf('Parallel child production runs: %d process workers\n',pool.NumWorkers);
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~)tick());
        parfor k=1:nTasks
            childBundles{k}=run_child(taskMode(k),taskCase(k,:),o,childRoot,k);
            send(Q,1);
        end
    else
        fprintf('Parallel child production runs: disabled\n');
        for k=1:nTasks
            childBundles{k}=run_child(taskMode(k),taskCase(k,:),o,childRoot,k);
            fprintf('  completed child %d/%d\n',k,nTasks);
        end
    end

    raw=build_raw_table(childBundles,taskMode,taskCase,o);
    [pairs,summary]=build_summary(raw,o);

    bundle=struct();
    bundle.meta=runMeta;
    bundle.experiment='i10_final_sa_pilot';
    bundle.run=struct('id',runId,'label',char(string(o.RunLabel)),'outputDirectory',outDir);
    bundle.cases=cases; bundle.boundaryCase=boundary;
    bundle.policy=struct('name','case-adaptive','ruleVersion','I9-v1', ...
        'baseDt',double(o.BaseDt),'refinedDt',double(o.RefinedDt), ...
        'selectedDt',selectedDt,'selectedMeta',{selectedMeta}, ...
        'referenceBoundaryDt',double(o.ReferenceBoundaryDt));
    bundle.method=struct('fastFadingMethod','logh','loghSpanSigma',double(o.loghSpanSigma), ...
        'loghYBlockSize',double(o.loghYBlockSize),'yGridMode','candidate-adaptive', ...
        'adaptiveScale',double(o.AdaptiveScale));
    bundle.sa=struct('maxIter',double(o.saMaxIter),'nStarts',double(o.saNStarts), ...
        'nReplicates',double(o.nReplicates),'T0',double(o.T0),'Tf',double(o.Tf), ...
        'baseStd0',double(o.BaseStd0),'itersPerTemp',double(o.ItersPerTemp), ...
        'baseSeed',double(o.baseSeed),'historyEvery',double(o.historyEvery));
    bundle.childRuns=cellfun(@child_descriptor,childBundles,'UniformOutput',false);
    bundle.results=raw; bundle.boundaryPairs=pairs; bundle.summary=summary;

    matPath=fullfile(outDir,'i10FinalSAPilot.mat'); save(matPath,'bundle','-v7.3');
    if o.WriteCSV
        writetable(raw,fullfile(outDir,'i10FinalSAPilot_raw.csv'));
        if ~isempty(pairs), writetable(pairs,fullfile(outDir,'i10FinalSAPilot_boundaryPairs.csv')); end
    end

    fprintf('\nSaved I.10 pilot bundle: %s\n',matPath);
    print_summary(summary,raw,pairs);

    function tick()
        doneCount=doneCount+1; fprintf('  completed child %d/%d\n',doneCount,nTasks);
    end
end


function b=run_child(mode,caseRow,o,childRoot,taskIndex)
    M=caseRow(1); SNR=caseRow(2); sig=caseRow(3);
    common={'MVec',M,'SNRVec',SNR,'TurbulenceVec',sig,'P_avg',o.P_avg, ...
        'minGap',o.minGap,'pinZero',o.pinZero,'FastFadingMethod','logh', ...
        'loghSpanSigma',o.loghSpanSigma,'loghYBlockSize',o.loghYBlockSize, ...
        'YGridMode','candidate-adaptive','YGridScale',o.AdaptiveScale, ...
        'saMaxIter',o.saMaxIter,'saNStarts',o.saNStarts,'nReplicates',o.nReplicates, ...
        'T0',o.T0,'Tf',o.Tf,'BaseStd0',o.BaseStd0,'ItersPerTemp',o.ItersPerTemp, ...
        'baseSeed',o.baseSeed,'historyEvery',o.historyEvery, ...
        'UseParallelCases',false,'MakePlots',false,'SaveFigures',false,'ResultsRoot',childRoot};

    if mode=="policy"
        runLabel=sprintf('I10_policy_%02d_M%d_SNR%s_sig%s',taskIndex,M,num_token(SNR,2),num_token(sig,4));
        b=sim_revision_AMI_vs_SNR(common{:},'LoghDtPolicy','case-adaptive', ...
            'loghDt',o.BaseDt,'loghDtRefined',o.RefinedDt,'RunLabel',runLabel);
    else
        runLabel=sprintf('I10_boundary_ref_M%d_SNR%s_sig%s',M,num_token(SNR,2),num_token(sig,4));
        b=sim_revision_AMI_vs_SNR(common{:},'LoghDtPolicy','fixed', ...
            'loghDt',o.ReferenceBoundaryDt,'loghDtRefined',o.RefinedDt,'RunLabel',runLabel);
    end
end


function T=build_raw_table(children,modes,cases,o)
    nTasks=numel(children); nRep=double(o.nReplicates); nRows=nTasks*nRep;
    Mode=strings(nRows,1); M=nan(nRows,1); SigmaX2=nan(nRows,1); SNRdB=nan(nRows,1);
    Replicate=nan(nRows,1); Seed=nan(nRows,1); SelectedDt=nan(nRows,1); IsRefined=false(nRows,1);
    AMIValidated=nan(nRows,1); AMIFast=nan(nRows,1); GainBits=nan(nRows,1); K99=nan(nRows,1);
    MaxAbsFastValidationGap=nan(nRows,1); SelectionChanged=false(nRows,1);
    RestartValidatedStd=nan(nRows,1); WallSeconds=nan(nRows,1); FastEvalPerSecond=nan(nRows,1);
    XMax=nan(nRows,1); WinnerX=cell(nRows,1); CaseFile=strings(nRows,1); ChildRunDir=strings(nRows,1);

    q=0;
    for k=1:nTasks
        b=children{k}; c=cases(k,:);
        for r=1:nRep
            q=q+1; Mode(q)=modes(k); M(q)=c(1); SNRdB(q)=c(2); SigmaX2(q)=c(3); Replicate(q)=r;
            Seed(q)=fso_result_utils.case_seed(double(o.baseSeed)+(r-1)*10000019,c(1),c(2),c(3));
            SelectedDt(q)=b.summary.raw.loghDtSelected(1,1,1,r);
            IsRefined(q)=logical(b.summary.raw.loghDtIsRefined(1,1,1,r));
            AMIValidated(q)=b.summary.raw.amiGSValidated(1,1,1,r);
            AMIFast(q)=b.summary.raw.amiGSFast(1,1,1,r);
            GainBits(q)=b.summary.raw.gainBits(1,1,1,r);
            K99(q)=b.summary.raw.convergenceIter99(1,1,1,r);
            MaxAbsFastValidationGap(q)=b.summary.raw.maxAbsFastValidationGap(1,1,1,r);
            SelectionChanged(q)=logical(b.summary.raw.selectionChanged(1,1,1,r));
            RestartValidatedStd(q)=b.summary.raw.restartValidatedStd(1,1,1,r);
            WallSeconds(q)=b.summary.raw.wallSeconds(1,1,1,r);
            FastEvalPerSecond(q)=b.summary.raw.fastEvaluationsPerSecond(1,1,1,r);
            XMax(q)=b.summary.raw.xMax(1,1,1,r);
            WinnerX{q}=b.summary.xWinner{1,1,1,r};
            CaseFile(q)=string(b.caseFiles{1,1,1,r});
            ChildRunDir(q)=string(b.run.outputDirectory);
        end
    end

    T=table(Mode,M,SigmaX2,SNRdB,Replicate,Seed,SelectedDt,IsRefined,AMIValidated,AMIFast, ...
        GainBits,K99,MaxAbsFastValidationGap,SelectionChanged,RestartValidatedStd,WallSeconds, ...
        FastEvalPerSecond,XMax,WinnerX,CaseFile,ChildRunDir);
end


function [P,S]=build_summary(T,o)
    policy=T(T.Mode=="policy",:);
    S=struct();
    S.nPolicyRows=height(policy); S.nReferenceRows=sum(T.Mode=="boundary-reference");
    S.fastValidationTolerance=double(o.FastValidationTolerance);
    S.boundaryValidatedDeltaTolerance=double(o.BoundaryValidatedDeltaTolerance);
    S.maxPolicyFastValidationGap=max(policy.MaxAbsFastValidationGap,[],'omitnan');
    S.meanPolicyFastValidationGap=mean(policy.MaxAbsFastValidationGap,'omitnan');
    S.nPolicyFastValidationFailures=sum(policy.MaxAbsFastValidationGap>o.FastValidationTolerance);
    S.nPolicySelectionChanges=sum(policy.SelectionChanged);
    S.selectionChangedFraction=mean(double(policy.SelectionChanged),'omitnan');
    S.maxPolicyK99=max(policy.K99,[],'omitnan');
    S.meanPolicyK99=mean(policy.K99,'omitnan');
    S.nPolicyK99AtBudget=sum(policy.K99>=double(o.saMaxIter));
    S.meanPolicyRestartValidatedStd=mean(policy.RestartValidatedStd,'omitnan');
    S.maxPolicyRestartValidatedStd=max(policy.RestartValidatedStd,[],'omitnan');
    S.totalPolicyWallSeconds=sum(policy.WallSeconds,'omitnan');
    S.meanPolicyFastEvalPerSecond=mean(policy.FastEvalPerSecond,'omitnan');

    P=table();
    if o.RunBoundaryReference
        boundary=double(o.BoundaryCase); nRep=double(o.nReplicates);
        Replicate=(1:nRep).'; Seed=nan(nRep,1); PolicyDt=nan(nRep,1); ReferenceDt=nan(nRep,1);
        PolicyAMI=nan(nRep,1); ReferenceAMI=nan(nRep,1); ValidatedDelta=nan(nRep,1); AbsValidatedDelta=nan(nRep,1);
        PolicyGain=nan(nRep,1); ReferenceGain=nan(nRep,1); PolicyK99=nan(nRep,1); ReferenceK99=nan(nRep,1);
        PolicyRestartStd=nan(nRep,1); ReferenceRestartStd=nan(nRep,1); PolicyWallSeconds=nan(nRep,1);
        ReferenceWallSeconds=nan(nRep,1); ReferenceOverPolicyWallRatio=nan(nRep,1);
        MaxAbsConstellationDelta=nan(nRep,1); SameSeed=false(nRep,1);
        for r=1:nRep
            key=T.M==boundary(1)&T.SNRdB==boundary(2)&abs(T.SigmaX2-boundary(3))<1e-12&T.Replicate==r;
            a=T(key&T.Mode=="policy",:); b=T(key&T.Mode=="boundary-reference",:);
            if height(a)~=1||height(b)~=1, error('sim_i10_final_sa_pilot:BoundaryPairing','Missing boundary A/B pair.'); end
            Seed(r)=a.Seed; SameSeed(r)=a.Seed==b.Seed; PolicyDt(r)=a.SelectedDt; ReferenceDt(r)=b.SelectedDt;
            PolicyAMI(r)=a.AMIValidated; ReferenceAMI(r)=b.AMIValidated;
            ValidatedDelta(r)=a.AMIValidated-b.AMIValidated; AbsValidatedDelta(r)=abs(ValidatedDelta(r));
            PolicyGain(r)=a.GainBits; ReferenceGain(r)=b.GainBits; PolicyK99(r)=a.K99; ReferenceK99(r)=b.K99;
            PolicyRestartStd(r)=a.RestartValidatedStd; ReferenceRestartStd(r)=b.RestartValidatedStd;
            PolicyWallSeconds(r)=a.WallSeconds; ReferenceWallSeconds(r)=b.WallSeconds;
            ReferenceOverPolicyWallRatio(r)=b.WallSeconds/max(a.WallSeconds,eps);
            xa=a.WinnerX{1}; xb=b.WinnerX{1};
            if numel(xa)==numel(xb), MaxAbsConstellationDelta(r)=max(abs(xa(:)-xb(:))); end
        end
        P=table(Replicate,Seed,PolicyDt,ReferenceDt,PolicyAMI,ReferenceAMI,ValidatedDelta,AbsValidatedDelta, ...
            PolicyGain,ReferenceGain,PolicyK99,ReferenceK99,PolicyRestartStd,ReferenceRestartStd, ...
            PolicyWallSeconds,ReferenceWallSeconds,ReferenceOverPolicyWallRatio,MaxAbsConstellationDelta,SameSeed);
        S.maxBoundaryAbsValidatedDelta=max(P.AbsValidatedDelta,[],'omitnan');
        S.meanBoundaryValidatedDelta=mean(P.ValidatedDelta,'omitnan');
        S.maxBoundaryConstellationDelta=max(P.MaxAbsConstellationDelta,[],'omitnan');
        S.meanBoundaryReferenceOverPolicyWallRatio=mean(P.ReferenceOverPolicyWallRatio,'omitnan');
        S.allBoundarySeedsMatch=all(P.SameSeed);
        S.nBoundaryValidatedDeltaFailures=sum(P.AbsValidatedDelta>o.BoundaryValidatedDeltaTolerance);
    else
        S.maxBoundaryAbsValidatedDelta=NaN; S.meanBoundaryValidatedDelta=NaN;
        S.maxBoundaryConstellationDelta=NaN; S.meanBoundaryReferenceOverPolicyWallRatio=NaN;
        S.allBoundarySeedsMatch=true; S.nBoundaryValidatedDeltaFailures=0;
    end

    S.policyFastValidationPass=S.nPolicyFastValidationFailures==0;
    S.policySelectionPass=S.nPolicySelectionChanges==0;
    S.boundaryReferencePass=S.nBoundaryValidatedDeltaFailures==0 && S.allBoundarySeedsMatch;
    S.overallPass=S.policyFastValidationPass && S.policySelectionPass && S.boundaryReferencePass;
end


function d=child_descriptor(b)
    d=struct('runId',b.run.id,'outputDirectory',b.run.outputDirectory, ...
        'loghDtPolicy',b.method.loghDtPolicy,'loghDtGrid',b.method.loghDtGrid, ...
        'nRefinedPhysicalCases',b.method.nRefinedPhysicalCases);
end


function print_summary(S,T,P)
    fprintf('\nI.10 FINAL SA PILOT SUMMARY\n');
    fprintf('Policy replicate rows                  : %d\n',S.nPolicyRows);
    fprintf('Reference replicate rows               : %d\n',S.nReferenceRows);
    fprintf('Policy max/mean F-V restart gap        : %.3e / %.3e\n',S.maxPolicyFastValidationGap,S.meanPolicyFastValidationGap);
    fprintf('Policy F-V failures (tol %.1e)         : %d\n',S.fastValidationTolerance,S.nPolicyFastValidationFailures);
    fprintf('Policy selection changes               : %d\n',S.nPolicySelectionChanges);
    fprintf('Policy k99 mean/max                    : %.1f / %.0f of budget\n',S.meanPolicyK99,S.maxPolicyK99);
    fprintf('Policy k99 at full budget              : %d\n',S.nPolicyK99AtBudget);
    fprintf('Policy restart validated std mean/max  : %.3e / %.3e\n',S.meanPolicyRestartValidatedStd,S.maxPolicyRestartValidatedStd);
    fprintf('Total policy wall time                 : %.2f h\n',S.totalPolicyWallSeconds/3600);
    if ~isempty(P)
        fprintf('Boundary max |policy-ref validated AMI|: %.3e (tol %.1e)\n', ...
            S.maxBoundaryAbsValidatedDelta,S.boundaryValidatedDeltaTolerance);
        fprintf('Boundary same-seed pairs               : %d\n',S.allBoundarySeedsMatch);
        fprintf('Boundary mean ref/policy wall ratio    : %.3fx\n',S.meanBoundaryReferenceOverPolicyWallRatio);
    end
    fprintf('OVERALL PILOT PASS                     : %d\n',S.overallPass);

    fprintf('\nPolicy rows:\n');
    disp(T(T.Mode=="policy",{'M','SigmaX2','SNRdB','Replicate','SelectedDt','AMIValidated', ...
        'K99','MaxAbsFastValidationGap','SelectionChanged','RestartValidatedStd','WallSeconds'}));
    if ~isempty(P)
        fprintf('\nBoundary paired comparison:\n'); disp(P);
    end
end


function tf=cases_valid(v)
    tf=isnumeric(v)&&ismatrix(v)&&size(v,2)==3&&~isempty(v)&&all(isfinite(v(:))) ...
        &&all(v(:,1)>=2)&&all(mod(v(:,1),1)==0)&&all(v(:,3)>=0);
end
function t=ternary(c,a,b), if c,t=a;else,t=b;end,end
function t=sanitize_label(v)
    t=strtrim(char(string(v))); if isempty(t),t='run';return;end
    t=regexprep(t,'[^A-Za-z0-9._-]+','_'); t=regexprep(t,'_+','_');
    t=regexprep(t,'^[._-]+|[._-]+$',''); if isempty(t),t='run';end
end
function [d,id]=make_unique_dir(parent,id)
    if ~exist(parent,'dir'),mkdir(parent);end
    d=fullfile(parent,id); base=id; n=1;
    while exist(d,'dir'),n=n+1;id=sprintf('%s_%02d',base,n);d=fullfile(parent,id);end
    [ok,msg]=mkdir(d); if ~ok,error('sim_i10_final_sa_pilot:mkdir','%s',msg);end
end
function token=num_token(v,nDec)
    token=sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token=strrep(token,'-','m'); token=strrep(token,'+','p'); token=strrep(token,'.','p');
end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=nonnegative_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;end
function tf=positive_integer(v),tf=positive_scalar(v)&&mod(v,1)==0;end
