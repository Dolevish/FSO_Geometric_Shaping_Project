function result = sim_k_min_gap_sensitivity(varargin)
%SIM_K_MIN_GAP_SENSITIVITY  Canonical Commit-K d_min sensitivity rerun.
%
% Frozen study:
%   M             = [16 32]
%   SNR           = 20 dB
%   sigma_R^2     = 0.2
%   d_min         = [0 0.005 0.01 0.02 0.05]
%   replicates    = 2
%   restarts      = 8
%   SA iterations = 8000/restart for BOTH M values
%
% All other numerical and SA settings are inherited from the final I.11.1
% production configuration, including the accurate log-h fast evaluator,
% candidate-adaptive y grid and independent AMI validator.
%
% A dynamic parfeval queue keeps one validated SA restart per process worker
% and immediately dispatches the next pending restart on completion.
%
% New results are written ONLY beneath results/ieee_revision2.
%
% New run:
%   Rk = sim_k_min_gap_sensitivity();
%
% Dry-run:
%   Rk = sim_k_min_gap_sensitivity('DryRun',true);
%
% Resume:
%   Rk = sim_k_min_gap_sensitivity('ResumeDirectory','/path/to/K_...');
%
% No global-optimum claim is implied.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ResultsRoot',fso_result_utils.revision2_results_root(),@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'ResumeDirectory','',@text_scalar);
    addParameter(p,'ParallelWorkers',8,@positive_integer);
    addParameter(p,'FreshPool',true,@logical_scalar);
    addParameter(p,'DeletePoolOnFinish',true,@logical_scalar);
    addParameter(p,'WriteCSV',true,@logical_scalar);
    addParameter(p,'MakeFigure',true,@logical_scalar);
    addParameter(p,'DryRun',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    K=k_min_gap_config();
    K.ParallelWorkers=double(o.ParallelWorkers);

    [casePlan,taskPlan]=k_build_task_plan(K);
    signature=build_signature(K,casePlan);
    signatureText=char(jsonencode(signature));

    fprintf('\n============================================================\n');
    fprintf('COMMIT K - MINIMUM-SPACING SENSITIVITY (%s)\n',K.version);
    fprintf('M=%s | SNR=%g dB | sigma_R^2=%.1f\n', ...
        mat2str(K.MVec),K.SNRdB,K.SigmaR2);
    fprintf('d_min=%s | replicates=%d | starts=%d | iter/start=%d\n', ...
        mat2str(K.MinGapVec),numel(K.Replicates),K.NStarts,K.MaxIter);
    fprintf('Sensitivity points=%d | replicate cases=%d | restart tasks=%d\n', ...
        K.nPhysicalSensitivityPoints,K.nReplicateCases,K.nRestartTasks);
    fprintf('Total SA iterations=%d | nominal fast objective calls=%d\n', ...
        K.totalSAIterations,K.nominalFastObjectiveEvaluations);
    fprintf('Fast evaluator=%s | y-grid=%s | dt=%s\n', ...
        K.FastFadingMethod,K.YGridMode,mat2str(K.SelectedDtByM));
    fprintf('Scheduler=%s | workers=%d | nested parallelism=OFF\n', ...
        K.scheduler,K.ParallelWorkers);
    fprintf('Results root: %s\n',char(o.ResultsRoot));
    fprintf('No global-optimum claim is implied.\n');
    fprintf('============================================================\n');

    if o.DryRun
        result=struct('version',K.version,'dryRun',true,'config',K, ...
            'casePlan',casePlan,'taskPlan',taskPlan,'signature',signature);
        return;
    end

    [outDir,manifestPath,manifest,resumeRequested]=prepare_run_directory( ...
        o,K,signatureText);
    checkpointDir=fullfile(outDir,'checkpoints');
    if exist(checkpointDir,'dir')~=7,mkdir(checkpointDir);end
    fprintf('Run directory: %s\n',outDir);

    baselines=build_baselines(K,casePlan);
    caseSpecs=build_case_specs(K,casePlan,baselines);

    nTasks=height(taskPlan);
    taskResults=cell(nTasks,1);
    reused=false(nTasks,1);
    checkpointPaths=strings(nTasks,1);

    for j=1:nTasks
        cp=fullfile(checkpointDir,[char(taskPlan.IdentityKey(j)) '.mat']);
        checkpointPaths(j)=string(cp);
        identity=task_identity(taskPlan(j,:),casePlan);
        [ok,r]=load_valid_checkpoint(cp,signatureText,identity);
        if ok
            r.checkpointReused=true;
            taskResults{j}=r;
            reused(j)=true;
        end
    end

    pending=find(~reused);
    fprintf('Checkpoint reuse: %d/%d complete | pending=%d\n', ...
        sum(reused),nTasks,numel(pending));

    if isempty(pending)
        segment=empty_segment(K.ParallelWorkers);
        trace=table();
        fprintf('All restart checkpoints already exist; skipping parallel execution.\n');
    else
        [taskResults,segment,trace]=run_dynamic_queue( ...
            taskResults,pending,taskPlan,casePlan,caseSpecs,checkpointPaths, ...
            signatureText,K,o);
        manifest.schedulerSegments(end+1)=segment;
        manifest.lastUpdatedUTC=utc_now();
        save_manifest_atomic(manifestPath,manifest);
    end

    if any(cellfun(@isempty,taskResults))
        error('sim_k_min_gap_sensitivity:Incomplete','One or more restart tasks are missing.');
    end

    [restartSummary,replicateSummary,physicalSummary,metrics]= ...
        reduce_results(taskResults,taskPlan,casePlan,baselines,K);

    metrics.checkpointReused=sum(reused);
    metrics.checkpointExecuted=nTasks-sum(reused);
    metrics.resumeRequested=resumeRequested;
    metrics.scheduler=aggregate_scheduler_metrics(manifest.schedulerSegments,K.ParallelWorkers);

    result=struct();
    result.version=K.version;
    result.config=K;
    result.signature=signature;
    result.casePlan=casePlan;
    result.taskPlan=taskPlan;
    result.baselines=baselines;
    result.restartSummary=restartSummary;
    result.replicateSummary=replicateSummary;
    result.physicalSummary=physicalSummary;
    result.metrics=metrics;
    result.schedulerTrace=trace;
    result.outputDirectory=outDir;
    result.detailedTaskResults=taskResults;
    result.globalOptimumClaim=false;

    matPath=fullfile(outDir,'K_min_gap_results.mat');
    save(matPath,'result','-v7.3');

    if o.WriteCSV
        writetable(restartSummary,fullfile(outDir,'K_min_gap_restart_results.csv'));
        writetable(replicateSummary,fullfile(outDir,'K_min_gap_replicate_summary.csv'));
        writetable(physicalSummary,fullfile(outDir,'K_min_gap_physical_summary.csv'));
        if ~isempty(trace)
            writetable(trace,fullfile(outDir,'K_scheduler_trace.csv'));
        end
    end

    if o.MakeFigure
        result.figurePaths=plot_k_min_gap_sensitivity(result,'OutputDirectory',outDir);
        save(matPath,'result','-v7.3');
    else
        result.figurePaths=struct();
    end

    manifest.completed=true;
    manifest.completedUTC=utc_now();
    manifest.lastUpdatedUTC=manifest.completedUTC;
    manifest.resultPath=matPath;
    save_manifest_atomic(manifestPath,manifest);

    print_summary(result);
end


function baselines=build_baselines(K,casePlan)
    baselines=struct('M',cell(numel(K.MVec),1),'x',[],'amiFast',[], ...
        'amiValidated',[],'fastValidationGap',[],'selectedDt',[]);
    fprintf('Validating %d deterministic Uniform-PAM baselines...\n',numel(K.MVec));

    for iM=1:numel(K.MVec)
        M=K.MVec(iM);
        row=casePlan(find(casePlan.M==M,1),:);
        cfg=build_case_cfg(K,row,0);
        x=define_constellation(M,K.P_avg,true,"mean");
        AMI_functions.assert_constellation_feasible(x,cfg,1e-10);
        fast=cfg.AMI_Evaluator(x(:));
        val=cfg.AMI_Validator(x(:));

        baselines(iM).M=M;
        baselines(iM).x=x(:);
        baselines(iM).amiFast=fast;
        baselines(iM).amiValidated=val;
        baselines(iM).fastValidationGap=fast-val;
        baselines(iM).selectedDt=row.SelectedDt;
        fprintf('  M=%d | PAM fast=%.6f | validated=%.6f | gap=%+.3e\n', ...
            M,fast,val,fast-val);
    end
end


function specs=build_case_specs(K,casePlan,baselines)
    specs=cell(height(casePlan),1);
    for c=1:height(casePlan)
        R=casePlan(c,:);
        cfg=build_case_cfg(K,R,R.MinGap);
        b=baselines([baselines.M]==R.M);
        specs{c}=struct('cfg',cfg,'xBaseline',b.x(:), ...
            'seed',uint32(R.Seed),'baselineValidated',b.amiValidated);
    end
end


function cfg=build_case_cfg(K,R,minGap)
    xmax=fso_result_utils.feasible_xmax_bound(R.M,K.P_avg,minGap,K.pinZero);
    cfg=build_fso_config(R.M,K.P_avg,R.SNRdB,R.SigmaR2, ...
        'FastFadingMethod',K.FastFadingMethod, ...
        'loghDt',R.SelectedDt, ...
        'loghSpanSigma',K.loghSpanSigma, ...
        'loghYBlockSize',K.loghYBlockSize, ...
        'YGridMode',K.YGridMode,'YGridScale',K.YGridScale, ...
        'xMaxBound',xmax, ...
        'saMaxIter',K.MaxIter,'saNStarts',K.NStarts, ...
        'saUseParallel',false,'minGap',minGap,'pinZero',K.pinZero, ...
        'seedInit',R.Seed,'logEvery',0,'historyEvery',K.historyEvery);

    % Explicitly re-freeze the production SA profile.
    cfg.SA.T0=K.T0;
    cfg.SA.Tf=K.Tf;
    cfg.SA.baseStd0=K.BaseStd0;
    cfg.SA.itersPerTemp=K.ItersPerTemp;
    cfg.SA.nBlocks=ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate=exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);
end


function sig=build_signature(K,casePlan)
    sig=struct();
    sig.version=K.version;
    sig.MVec=K.MVec;
    sig.SNRdB=K.SNRdB;
    sig.SigmaR2=K.SigmaR2;
    sig.MinGapVec=K.MinGapVec;
    sig.Replicates=K.Replicates;
    sig.NStarts=K.NStarts;
    sig.MaxIter=K.MaxIter;
    sig.P_avg=K.P_avg;
    sig.pinZero=K.pinZero;
    sig.FastFadingMethod=K.FastFadingMethod;
    sig.LoghDtPolicy=K.LoghDtPolicy;
    sig.SelectedDtByM=K.SelectedDtByM;
    sig.loghSpanSigma=K.loghSpanSigma;
    sig.loghYBlockSize=K.loghYBlockSize;
    sig.YGridMode=K.YGridMode;
    sig.YGridScale=K.YGridScale;
    sig.T0=K.T0;
    sig.Tf=K.Tf;
    sig.BaseStd0=K.BaseStd0;
    sig.ItersPerTemp=K.ItersPerTemp;
    sig.historyEvery=K.historyEvery;
    sig.baseSeed=K.baseSeed;
    sig.scheduler=K.scheduler;
    sig.caseIdentities=table2struct(casePlan(:, ...
        {'M','SNRdB','SigmaR2','MinGap','Replicate','Seed','NStarts','MaxIter','SelectedDt'}));
end


function [outDir,manifestPath,manifest,resumeRequested]=prepare_run_directory(o,K,signatureText)
    resumeRequested=strlength(string(o.ResumeDirectory))>0;

    if resumeRequested
        outDir=char(string(o.ResumeDirectory));
        manifestPath=fullfile(outDir,'K_run_manifest.mat');
        if exist(manifestPath,'file')~=2
            error('sim_k_min_gap_sensitivity:ManifestMissing', ...
                'Resume directory has no K_run_manifest.mat: %s',outDir);
        end
        X=load(manifestPath,'manifest');
        manifest=X.manifest;
        if ~isfield(manifest,'signatureText')||~strcmp(manifest.signatureText,signatureText)
            error('sim_k_min_gap_sensitivity:ResumeMismatch', ...
                'Resume directory does not match the frozen Commit-K signature.');
        end
        fprintf('Resuming Commit-K run: %s\n',outDir);
        return;
    end

    if strlength(string(o.OutputDirectory))>0
        parent=char(string(o.OutputDirectory));
    else
        parent=fullfile(char(o.ResultsRoot),'commit_k_min_gap');
    end
    ensure_under_revision2(parent);

    if exist(parent,'dir')~=7,mkdir(parent);end
    token=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    outDir=fullfile(parent,['K_' token]);
    [ok,msg]=mkdir(outDir);
    if ~ok,error('sim_k_min_gap_sensitivity:mkdir','%s',msg);end

    manifestPath=fullfile(outDir,'K_run_manifest.mat');
    manifest=struct();
    manifest.version=K.version;
    manifest.signatureText=signatureText;
    manifest.createdUTC=utc_now();
    manifest.lastUpdatedUTC=manifest.createdUTC;
    manifest.completed=false;
    manifest.schedulerSegments=repmat(empty_segment(K.ParallelWorkers),0,1);
    manifest.resultPath='';
    save_manifest_atomic(manifestPath,manifest);
end


function ensure_under_revision2(path)
    canonical=fso_result_utils.revision2_results_root();
    p=char(java.io.File(path).getCanonicalPath());
    r=char(java.io.File(canonical).getCanonicalPath());
    if ~(strcmp(p,r)||startsWith(p,[r filesep]))
        error('sim_k_min_gap_sensitivity:ResultsRootPolicy', ...
            ['Commit K results must be written under results/ieee_revision2. ' ...
             'Requested path: %s'],path);
    end
end


function [taskResults,segment,trace]=run_dynamic_queue( ...
    taskResults,pending,taskPlan,casePlan,caseSpecs,checkpointPaths, ...
    signatureText,K,o)

    if license('test','Distrib_Computing_Toolbox')==0
        error('sim_k_min_gap_sensitivity:ParallelToolbox', ...
            'Parallel Computing Toolbox is required for the Commit-K final run.');
    end

    pp=gcp('nocreate');
    if ~isempty(pp)&&o.FreshPool
        delete(pp); pp=[];
    end
    if isempty(pp)
        pp=parpool('Processes',double(o.ParallelWorkers));
    elseif pp.NumWorkers~=double(o.ParallelWorkers)
        delete(pp);
        pp=parpool('Processes',double(o.ParallelWorkers));
    end
    if pp.NumWorkers~=double(o.ParallelWorkers)
        error('sim_k_min_gap_sensitivity:PoolSize', ...
            'Requested %d workers, got %d.',o.ParallelWorkers,pp.NumWorkers);
    end
    poolCleanup=onCleanup(@() cleanup_pool(pp,o.DeletePoolOnFinish)); %#ok<NASGU>

    nWorkers=pp.NumWorkers;
    nPending=numel(pending);
    active=parallel.FevalFuture.empty(0,1);
    activeRows=zeros(0,1);
    nextPending=1;
    completed=0;

    traceElapsed=zeros(nPending,1);
    traceCompleted=zeros(nPending,1);
    traceActive=zeros(nPending,1);
    tracePending=zeros(nPending,1);
    traceIdentity=strings(nPending,1);
    traceWorkerTaskSeconds=zeros(nPending,1);

    wallTimer=tic;

    while numel(active)<nWorkers && nextPending<=nPending
        rowIdx=pending(nextPending);
        active(end+1,1)=submit_one(pp,rowIdx,taskPlan,casePlan,caseSpecs, ...
            checkpointPaths,signatureText); %#ok<AGROW>
        activeRows(end+1,1)=rowIdx; %#ok<AGROW>
        nextPending=nextPending+1;
    end

    try
        while ~isempty(active)
            [finishedIdx,res]=fetchNext(active);
            rowIdx=activeRows(finishedIdx);
            active(finishedIdx)=[];
            activeRows(finishedIdx)=[];

            taskResults{rowIdx}=res;
            completed=completed+1;

            if nextPending<=nPending
                newRow=pending(nextPending);
                active(end+1,1)=submit_one(pp,newRow,taskPlan,casePlan,caseSpecs, ...
                    checkpointPaths,signatureText);
                activeRows(end+1,1)=newRow;
                nextPending=nextPending+1;
            end

            traceElapsed(completed)=toc(wallTimer);
            traceCompleted(completed)=completed;
            traceActive(completed)=numel(active);
            tracePending(completed)=max(0,nPending-nextPending+1);
            traceIdentity(completed)=taskPlan.IdentityKey(rowIdx);
            traceWorkerTaskSeconds(completed)=res.taskWallSeconds;

            fprintf('  K completed %d/%d this run | total %d/%d | active=%d | pending=%d | %s\n', ...
                completed,nPending,sum(~cellfun(@isempty,taskResults)),height(taskPlan), ...
                numel(active),tracePending(completed),taskPlan.IdentityKey(rowIdx));
        end
    catch ME
        if ~isempty(active),try,cancel(active);catch,end,end
        rethrow(ME);
    end

    wallSeconds=toc(wallTimer);
    busySeconds=sum(traceWorkerTaskSeconds(1:completed));
    capacitySeconds=nWorkers*wallSeconds;
    utilization=busySeconds/max(capacitySeconds,eps);

    trace=table(traceElapsed(1:completed),traceCompleted(1:completed), ...
        traceActive(1:completed),tracePending(1:completed), ...
        traceIdentity(1:completed),traceWorkerTaskSeconds(1:completed), ...
        'VariableNames',{'ElapsedSeconds','CompletedThisRun','ActiveAfterDispatch', ...
        'PendingAfterDispatch','FinishedIdentity','FinishedTaskWallSeconds'});

    preTail=trace.PendingAfterDispatch>0;
    if any(preTail)
        fullQueueMaintained=all(trace.ActiveAfterDispatch(preTail)==nWorkers);
        minActiveBeforeTail=min(trace.ActiveAfterDispatch(preTail));
    else
        fullQueueMaintained=true;
        minActiveBeforeTail=nWorkers;
    end

    segment=struct();
    segment.completedUTC=utc_now();
    segment.workers=nWorkers;
    segment.executedTasks=completed;
    segment.wallSeconds=wallSeconds;
    segment.workerBusySeconds=busySeconds;
    segment.workerCapacitySeconds=capacitySeconds;
    segment.workerUtilization=utilization;
    segment.fullQueueMaintainedUntilTail=fullQueueMaintained;
    segment.minActiveBeforeTail=minActiveBeforeTail;
    segment.scheduler=K.scheduler;

    fprintf('\nCommit-K scheduler: wall=%.2f h | busy=%.2f worker-h | utilization=%.1f%%\n', ...
        wallSeconds/3600,busySeconds/3600,100*utilization);
    fprintf('Dynamic queue full before tail: %d | min active while pending>0: %d/%d\n', ...
        fullQueueMaintained,minActiveBeforeTail,nWorkers);
end


function fut=submit_one(pp,rowIdx,taskPlan,casePlan,caseSpecs,checkpointPaths,signatureText)
    t=taskPlan(rowIdx,:);
    c=t.TaskCaseIndex;
    R=casePlan(c,:);
    spec=caseSpecs{c};
    identity=task_identity(t,casePlan);
    label=sprintf('[K|M%d|d%.3f|rep%d]',R.M,R.MinGap,R.Replicate);

    fut=parfeval(pp,@sa_restart_task_checkpointed,1, ...
        spec.cfg,spec.xBaseline,spec.seed,t.Restart,label, ...
        char(checkpointPaths(rowIdx)),signatureText,identity);
end


function id=task_identity(t,casePlan)
    R=casePlan(t.TaskCaseIndex,:);
    id=struct('IdentityKey',char(t.IdentityKey),'M',R.M,'SNRdB',R.SNRdB, ...
        'SigmaR2',R.SigmaR2,'MinGap',R.MinGap,'Replicate',R.Replicate, ...
        'Restart',t.Restart,'Seed',double(R.Seed));
end


function [ok,r]=load_valid_checkpoint(path,signatureText,identity)
    ok=false; r=[];
    if exist(path,'file')~=2,return;end
    try
        S=load(path,'restartResult','checkpointMeta');
        ok=isfield(S,'restartResult')&&isfield(S,'checkpointMeta')&& ...
            isfield(S.checkpointMeta,'signatureText')&& ...
            isfield(S.checkpointMeta,'identity')&& ...
            strcmp(S.checkpointMeta.signatureText,signatureText)&& ...
            isequaln(S.checkpointMeta.identity,identity);
        if ok,r=S.restartResult;end
    catch
        ok=false; r=[];
    end
end


function [restartSummary,replicateSummary,physicalSummary,metrics]= ...
        reduce_results(taskResults,taskPlan,casePlan,baselines,K)

    n=height(taskPlan);
    M=zeros(n,1); SNRdB=zeros(n,1); SigmaR2=zeros(n,1); MinGap=zeros(n,1);
    Replicate=zeros(n,1); Restart=zeros(n,1); Seed=zeros(n,1);
    FastAMI=zeros(n,1); ValidatedAMI=zeros(n,1); FastValidatorGap=zeros(n,1);
    InitialFastAMI=zeros(n,1); AcceptanceRate=zeros(n,1);
    SARuntimeSeconds=zeros(n,1); ValidationSeconds=zeros(n,1);
    TaskWallSeconds=zeros(n,1); ActualMinGap=zeros(n,1); XMax=zeros(n,1);
    IdentityKey=strings(n,1); CheckpointReused=false(n,1);

    for j=1:n
        t=taskPlan(j,:);
        R=casePlan(t.TaskCaseIndex,:);
        z=taskResults{j};

        M(j)=R.M; SNRdB(j)=R.SNRdB; SigmaR2(j)=R.SigmaR2; MinGap(j)=R.MinGap;
        Replicate(j)=R.Replicate; Restart(j)=t.Restart; Seed(j)=R.Seed;
        FastAMI(j)=z.bestMIFast; ValidatedAMI(j)=z.bestMIValidated;
        FastValidatorGap(j)=z.validationGap; InitialFastAMI(j)=z.initialMI;
        AcceptanceRate(j)=z.acceptanceRate; SARuntimeSeconds(j)=z.runtime;
        ValidationSeconds(j)=z.validationRuntime; TaskWallSeconds(j)=z.taskWallSeconds;
        ActualMinGap(j)=safe_min_diff(z.x_best); XMax(j)=max(z.x_best);
        IdentityKey(j)=t.IdentityKey;
        if isfield(z,'checkpointReused'),CheckpointReused(j)=logical(z.checkpointReused);end
    end

    restartSummary=table(M,SNRdB,SigmaR2,MinGap,Replicate,Restart,Seed,IdentityKey, ...
        InitialFastAMI,FastAMI,ValidatedAMI,FastValidatorGap,AcceptanceRate, ...
        SARuntimeSeconds,ValidationSeconds,TaskWallSeconds,ActualMinGap,XMax,CheckpointReused);
    restartSummary=sortrows(restartSummary,{'M','MinGap','Replicate','Restart'});

    nc=height(casePlan);
    WinnerRestart=zeros(nc,1); WinnerFastAMI=zeros(nc,1); WinnerValidatedAMI=zeros(nc,1);
    PAMValidatedAMI=zeros(nc,1); GainOverPAM=zeros(nc,1);
    RestartValidatedMean=zeros(nc,1); RestartValidatedStd=zeros(nc,1);
    RestartValidatedMin=zeros(nc,1); RestartValidatedMax=zeros(nc,1);
    RestartValidatedRange=zeros(nc,1); MaxAbsFastValidatorGap=zeros(nc,1);
    FastWinnerRestart=zeros(nc,1); SelectionChanged=false(nc,1);
    WinnerActualMinGap=zeros(nc,1); WinnerXMax=zeros(nc,1); WinnerConstellation=cell(nc,1);

    for c=1:nc
        C=casePlan(c,:);
        idx=restartSummary.M==C.M & abs(restartSummary.MinGap-C.MinGap)<1e-12 & ...
            restartSummary.Replicate==C.Replicate;
        R=restartSummary(idx,:);
        if height(R)~=K.NStarts
            error('sim_k_min_gap_sensitivity:RestartReduction', ...
                'Expected %d restarts for case %d, got %d.',K.NStarts,c,height(R));
        end

        [WinnerValidatedAMI(c),iv]=max(R.ValidatedAMI);
        [~,ifast]=max(R.FastAMI);
        WinnerRestart(c)=R.Restart(iv);
        FastWinnerRestart(c)=R.Restart(ifast);
        SelectionChanged(c)=WinnerRestart(c)~=FastWinnerRestart(c);
        WinnerFastAMI(c)=R.FastAMI(iv);

        b=baselines([baselines.M]==C.M);
        PAMValidatedAMI(c)=b.amiValidated;
        GainOverPAM(c)=WinnerValidatedAMI(c)-PAMValidatedAMI(c);

        RestartValidatedMean(c)=mean(R.ValidatedAMI);
        RestartValidatedStd(c)=std(R.ValidatedAMI,0);
        RestartValidatedMin(c)=min(R.ValidatedAMI);
        RestartValidatedMax(c)=max(R.ValidatedAMI);
        RestartValidatedRange(c)=RestartValidatedMax(c)-RestartValidatedMin(c);
        MaxAbsFastValidatorGap(c)=max(abs(R.FastValidatorGap));
        WinnerActualMinGap(c)=R.ActualMinGap(iv);
        WinnerXMax(c)=R.XMax(iv);

        rowInTask=find(taskPlan.TaskCaseIndex==c & taskPlan.Restart==WinnerRestart(c),1);
        WinnerConstellation{c}=taskResults{rowInTask}.x_best(:);
    end

    replicateSummary=[casePlan table(WinnerRestart,FastWinnerRestart,SelectionChanged, ...
        WinnerFastAMI,WinnerValidatedAMI,PAMValidatedAMI,GainOverPAM, ...
        RestartValidatedMean,RestartValidatedStd,RestartValidatedMin, ...
        RestartValidatedMax,RestartValidatedRange,MaxAbsFastValidatorGap, ...
        WinnerActualMinGap,WinnerXMax,WinnerConstellation)];
    replicateSummary=sortrows(replicateSummary,{'M','MinGap','Replicate'});

    keys=unique(replicateSummary(:,{'M','SNRdB','SigmaR2','MinGap'}),'rows','sorted');
    np=height(keys);
    Rep1ValidatedAMI=zeros(np,1); Rep2ValidatedAMI=zeros(np,1);
    MeanValidatedAMI=zeros(np,1); StdAcrossReplicates=zeros(np,1);
    ReplicateRange=zeros(np,1); PAMValidated=zeros(np,1); MeanGainOverPAM=zeros(np,1);
    GainRetentionVsGap0=zeros(np,1);

    for k=1:np
        idx=replicateSummary.M==keys.M(k) & ...
            abs(replicateSummary.MinGap-keys.MinGap(k))<1e-12;
        R=sortrows(replicateSummary(idx,:),'Replicate');
        if height(R)~=2||~isequal(R.Replicate(:).',[1 2])
            error('sim_k_min_gap_sensitivity:ReplicateReduction', ...
                'Expected exactly replicates 1 and 2.');
        end
        Rep1ValidatedAMI(k)=R.WinnerValidatedAMI(1);
        Rep2ValidatedAMI(k)=R.WinnerValidatedAMI(2);
        MeanValidatedAMI(k)=mean([Rep1ValidatedAMI(k) Rep2ValidatedAMI(k)]);
        StdAcrossReplicates(k)=std([Rep1ValidatedAMI(k) Rep2ValidatedAMI(k)],0);
        ReplicateRange(k)=abs(Rep1ValidatedAMI(k)-Rep2ValidatedAMI(k));
        PAMValidated(k)=mean(R.PAMValidatedAMI);
        MeanGainOverPAM(k)=MeanValidatedAMI(k)-PAMValidated(k);
    end

    physicalSummary=[keys table(PAMValidated,Rep1ValidatedAMI,Rep2ValidatedAMI, ...
        MeanValidatedAMI,StdAcrossReplicates,ReplicateRange,MeanGainOverPAM, ...
        GainRetentionVsGap0)];

    for m=K.MVec
        idxM=physicalSummary.M==m;
        idx0=idxM & abs(physicalSummary.MinGap)<1e-12;
        g0=physicalSummary.MeanGainOverPAM(idx0);
        if isempty(g0)||abs(g0)<eps
            physicalSummary.GainRetentionVsGap0(idxM)=NaN;
        else
            physicalSummary.GainRetentionVsGap0(idxM)= ...
                physicalSummary.MeanGainOverPAM(idxM)/g0;
        end
    end

    metrics=struct();
    metrics.nRestartTasks=height(restartSummary);
    metrics.nReplicateCases=height(replicateSummary);
    metrics.nPhysicalSensitivityPoints=height(physicalSummary);
    metrics.maxAbsFastValidatorGap=max(abs(restartSummary.FastValidatorGap));
    metrics.selectionChanges=sum(replicateSummary.SelectionChanged);
    metrics.maxReplicateRange=max(physicalSummary.ReplicateRange);
    metrics.meanAcceptanceRate=mean(restartSummary.AcceptanceRate);
    metrics.totalSAIterations=K.totalSAIterations;
    metrics.nominalFastObjectiveEvaluations=K.nominalFastObjectiveEvaluations;
end


function A=aggregate_scheduler_metrics(segments,workers)
    A=struct('nSegments',numel(segments),'wallSeconds',0,'workerBusySeconds',0, ...
        'workerCapacitySeconds',0,'workerUtilization',NaN, ...
        'allSegmentsMaintainedFullQueue',true,'minActiveBeforeTail',workers);
    if isempty(segments),return;end
    for k=1:numel(segments)
        A.wallSeconds=A.wallSeconds+segments(k).wallSeconds;
        A.workerBusySeconds=A.workerBusySeconds+segments(k).workerBusySeconds;
        A.workerCapacitySeconds=A.workerCapacitySeconds+segments(k).workerCapacitySeconds;
        A.allSegmentsMaintainedFullQueue=A.allSegmentsMaintainedFullQueue && ...
            segments(k).fullQueueMaintainedUntilTail;
        A.minActiveBeforeTail=min(A.minActiveBeforeTail,segments(k).minActiveBeforeTail);
    end
    A.workerUtilization=A.workerBusySeconds/max(A.workerCapacitySeconds,eps);
end


function print_summary(R)
    fprintf('\n============================================================\n');
    fprintf('COMMIT K COMPLETE\n');
    fprintf('Physical sensitivity points: %d | replicate cases: %d | restart tasks: %d\n', ...
        R.metrics.nPhysicalSensitivityPoints,R.metrics.nReplicateCases,R.metrics.nRestartTasks);
    fprintf('Max |fast-validator| gap: %.3e bit/symbol\n',R.metrics.maxAbsFastValidatorGap);
    fprintf('Fast-vs-validated winner changes: %d/%d replicate cases\n', ...
        R.metrics.selectionChanges,R.metrics.nReplicateCases);
    fprintf('Max replicate AMI range: %.6f bit/symbol\n',R.metrics.maxReplicateRange);
    if R.metrics.scheduler.nSegments>0
        fprintf('Worker utilization: %.1f%% | full queue before tail=%d | min active=%d/%d\n', ...
            100*R.metrics.scheduler.workerUtilization, ...
            R.metrics.scheduler.allSegmentsMaintainedFullQueue, ...
            R.metrics.scheduler.minActiveBeforeTail,R.config.ParallelWorkers);
    end
    fprintf('\nAMI vs d_min:\n');
    disp(R.physicalSummary(:,{'M','MinGap','Rep1ValidatedAMI','Rep2ValidatedAMI', ...
        'MeanValidatedAMI','StdAcrossReplicates','MeanGainOverPAM','GainRetentionVsGap0'}));
    fprintf('Output: %s\n',R.outputDirectory);
    fprintf('No global-optimum claim is made.\n');
    fprintf('============================================================\n');
end


function s=empty_segment(workers)
    s=struct('completedUTC','','workers',workers,'executedTasks',0, ...
        'wallSeconds',0,'workerBusySeconds',0,'workerCapacitySeconds',0, ...
        'workerUtilization',NaN,'fullQueueMaintainedUntilTail',true, ...
        'minActiveBeforeTail',workers,'scheduler','parfeval-dynamic-v1');
end

function save_manifest_atomic(path,manifest)
    folder=fileparts(path);
    tmp=[tempname(folder) '.mat'];
    cleanup=onCleanup(@() cleanup_tmp(tmp)); %#ok<NASGU>
    save(tmp,'manifest','-v7');
    [ok,msg]=movefile(tmp,path,'f');
    if ~ok,error('sim_k_min_gap_sensitivity:ManifestSave','%s',msg);end
end

function cleanup_tmp(path)
    if exist(path,'file')==2,try,delete(path);catch,end,end
end

function cleanup_pool(pp,doDelete)
    if doDelete&&~isempty(pp)
        try,if isvalid(pp),delete(pp);end,catch,end
    end
end

function d=safe_min_diff(x)
    x=sort(real(x(:)));
    if numel(x)<2,d=NaN;else,d=min(diff(x));end
end

function s=utc_now()
    s=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0&&mod(v,1)==0;end
