function result = sim_i13d_final_optimizer_benchmark(varargin)
%SIM_I13D_FINAL_OPTIMIZER_BENCHMARK
% Final Reviewer-2 SA-vs-Pattern-Search benchmark.
%
% This run is deliberately frozen before seeing the final results:
%   M          = [4 8 16 32]
%   SNR [dB]   = [20 30]
%   sigma_R^2  = [0 0.3]
%   replicates = [1 2]
%   starts     = frozen I.11.1 map [6 6 8 8]
%   PS budget  = SA iterations + 1 evaluations/start
%   tie        = |I_SA-I_PS| <= 1e-3 bit/symbol
%
% PARALLEL SCHEDULER
% ------------------
% This function intentionally does NOT use parfor.  It uses a dynamic
% parfeval work queue with restart-level tasks. At most one restart is in
% flight per worker. Whenever a worker finishes, fetchNext() receives that
% result and the next longest-first pending task is immediately submitted.
% Therefore, while at least one pending task remains, the scheduler keeps
% exactly Nworkers futures in flight (subject only to submit/fetch latency).
%
% There is no nested Pattern Search parallelism.
%
% CHECKPOINT / RESUME
% -------------------
% Each successful restart writes an atomic checkpoint on the worker.
% Resume validates a frozen scientific signature and reuses completed
% checkpoints.  Results are reduced only after every restart is available.
%
% New run:
%   R = sim_i13d_final_optimizer_benchmark();
%
% Resume:
%   R = sim_i13d_final_optimizer_benchmark('ResumeDirectory','/path/to/run');
%
% Dry-run plan:
%   R = sim_i13d_final_optimizer_benchmark('DryRun',true);
%
% No global-optimum claim is implied.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'ResumeDirectory','',@text_scalar);
    addParameter(p,'ParallelWorkers',8,@positive_integer);
    addParameter(p,'FreshPool',true,@logical_scalar);
    addParameter(p,'DeletePoolOnFinish',true,@logical_scalar);
    addParameter(p,'WriteCSV',true,@logical_scalar);
    addParameter(p,'MakeFigures',true,@logical_scalar);
    addParameter(p,'DryRun',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    if exist('patternsearch','file')~=2
        error('sim_i13d_final_optimizer_benchmark:ToolboxMissing', ...
            'Global Optimization Toolbox is required: patternsearch not found.');
    end
    if o.ParallelWorkers<1
        error('sim_i13d_final_optimizer_benchmark:Workers','ParallelWorkers must be positive.');
    end

    F=i13d_benchmark_config();
    F.ParallelWorkers=double(o.ParallelWorkers);
    F.FreshPool=logical(o.FreshPool);

    A=analyze_i11p1_production_results( ...
        'ProductionBundlePath',o.ProductionBundlePath, ...
        'ResultsRoot',o.ResultsRoot,'WriteCSV',true);
    S=load(A.sourceBundlePath,'bundle');
    if ~isfield(S,'bundle')
        error('sim_i13d_final_optimizer_benchmark:Bundle','Production bundle variable is missing.');
    end
    prod=S.bundle;

    [casePlan,taskPlan]=i13d_build_task_plan(A.caseTable,F);
    caseSpecs=build_case_specs(casePlan,prod,F);
    signature=build_signature(F,casePlan,A.sourceBundlePath);
    signatureText=char(jsonencode(signature));

    fprintf('\n============================================================\n');
    fprintf('I.13D FINAL OPTIMIZER BENCHMARK (%s)\n',F.version);
    fprintf('SA source: frozen %s production\n',F.sourceProductionFreeze);
    fprintf('Grid: M=%s | SNR=%s dB | sigma_R^2=%s | reps=%s\n', ...
        mat2str(F.MVec),mat2str(F.SNRVec),mat2str(F.TurbulenceVec),mat2str(F.Replicates));
    fprintf('Cases=%d | restart tasks=%d | starts=%s\n', ...
        F.nOptimizationCases,F.nRestartTasks,mat2str(F.RestartsByM));
    fprintf('Max PS objective evaluations across benchmark: %.0f\n',F.totalMaxFunctionEvaluations);
    fprintf('PS eval/start=%s | equivalence=%.1e bit/symbol\n', ...
        mat2str(F.MaxEvalsPerStart),F.tieToleranceBits);
    fprintf('Scheduler: %s | workers=%d | nested parallelism=OFF\n', ...
        F.scheduler,F.ParallelWorkers);
    fprintf('No global-optimum claim is implied.\n');
    fprintf('============================================================\n');

    if o.DryRun
        result=struct('version',F.version,'dryRun',true,'config',F, ...
            'sourceProduction',A.sourceBundlePath,'casePlan',casePlan,'taskPlan',taskPlan, ...
            'signature',signature);
        return;
    end

    [outDir,manifestPath,manifest,resumeRequested]=prepare_run_directory( ...
        o,A,F,signatureText);
    checkpointDir=fullfile(outDir,'checkpoints');
    if exist(checkpointDir,'dir')~=7,mkdir(checkpointDir);end

    % ------------------------------------------------------------------
    % Discover valid checkpoints before any pool is created.
    % taskResults is indexed by row in the priority-sorted taskPlan.
    % ------------------------------------------------------------------
    nTasks=height(taskPlan);
    taskResults=cell(nTasks,1);
    reused=false(nTasks,1);
    checkpointPaths=strings(nTasks,1);

    for j=1:nTasks
        cp=fullfile(checkpointDir,[char(taskPlan.IdentityKey(j)) '.mat']);
        checkpointPaths(j)=string(cp);
        if exist(cp,'file')==2
            C=load(cp,'checkpoint');
            validate_checkpoint(C.checkpoint,signatureText,taskPlan.IdentityKey(j));
            if ~C.checkpoint.result.success
                error('sim_i13d_final_optimizer_benchmark:BadCheckpoint', ...
                    'Checkpoint is not successful: %s',cp);
            end
            taskResults{j}=C.checkpoint.result;
            reused(j)=true;
        end
    end

    pending=find(~reused);
    fprintf('Checkpoint reuse: %d/%d complete | pending=%d\n',sum(reused),nTasks,numel(pending));

    if isempty(pending)
        fprintf('All restart checkpoints already exist; skipping parallel execution.\n');
        segment=empty_segment(F.ParallelWorkers);
        trace=table();
    else
        [taskResults,segment,trace]=run_dynamic_queue( ...
            taskResults,pending,taskPlan,casePlan,caseSpecs,checkpointPaths, ...
            signatureText,F,o);
        manifest.schedulerSegments(end+1)=segment;
        manifest.lastUpdatedUTC=utc_now();
        manifest.completed=false;
        save_manifest_atomic(manifestPath,manifest);
    end

    % Confirm all tasks are available and successful.
    for j=1:nTasks
        if isempty(taskResults{j}) || ~taskResults{j}.success
            error('sim_i13d_final_optimizer_benchmark:Incomplete', ...
                'Task %s is missing or unsuccessful.',taskPlan.IdentityKey(j));
        end
    end

    [restartSummary,caseSummary,physicalSummary,metrics]=reduce_results( ...
        taskResults,taskPlan,casePlan,F);

    % Optional Zhou-2025 baseline already generated in I.13B. It is joined
    % only if the saved result exists; the optimizer benchmark itself does
    % not depend on this file.
    [physicalSummary,zhouSource]=attach_zhou_if_available(physicalSummary,outDir);

    aggregate=aggregate_scheduler_metrics(manifest.schedulerSegments,F.ParallelWorkers);
    metrics.scheduler=aggregate;
    metrics.checkpointReused=sum(reused);
    metrics.checkpointExecuted=nTasks-sum(reused);
    metrics.resumeRequested=resumeRequested;

    result=struct();
    result.version=F.version;
    result.isReviewerFinal=true;
    result.sourceProduction=A.sourceBundlePath;
    result.config=F;
    result.signature=signature;
    result.casePlan=casePlan;
    result.taskPlan=taskPlan;
    result.restartSummary=restartSummary;
    result.caseSummary=caseSummary;
    result.physicalSummary=physicalSummary;
    result.metrics=metrics;
    result.schedulerTrace=trace;
    result.zhouSource=zhouSource;
    result.outputDirectory=outDir;
    result.detailedTaskResults=taskResults;
    result.globalOptimumClaim=false;

    matPath=fullfile(outDir,'i13d_optimizer_benchmark.mat');
    save(matPath,'result','-v7.3');

    if o.WriteCSV
        writetable(restartSummary,fullfile(outDir,'i13d_restart_tasks.csv'));
        writetable(caseSummary,fullfile(outDir,'i13d_case_summary.csv'));
        writetable(physicalSummary,fullfile(outDir,'i13d_physical_summary.csv'));
        if ~isempty(trace)
            writetable(trace,fullfile(outDir,'i13d_scheduler_trace.csv'));
        end
    end

    if o.MakeFigures
        result.figures=plot_i13d_optimizer_benchmark(result,'OutputDirectory',outDir);
        save(matPath,'result','-v7.3');
    else
        result.figures=struct();
    end

    manifest.completed=true;
    manifest.completedUTC=utc_now();
    manifest.resultPath=matPath;
    manifest.lastUpdatedUTC=manifest.completedUTC;
    save_manifest_atomic(manifestPath,manifest);

    print_summary(result);

end


function specs=build_case_specs(casePlan,prod,F)
    specs=cell(height(casePlan),1);
    P=i11p1_production_config();

    for c=1:height(casePlan)
        R=casePlan(c,:);
        [xBaseline,baselineAMI]=production_baseline(prod,R.M,R.SNRdB,R.SigmaR2);
        if abs(baselineAMI-R.PAMValidated)>1e-8
            error('sim_i13d_final_optimizer_benchmark:BaselineMismatch', ...
                'Production baseline mismatch for case %d.',c);
        end

        cfg=build_fso_config(R.M,F.P_avg,R.SNRdB,R.SigmaR2, ...
            'FastFadingMethod',P.FastFadingMethod, ...
            'loghDt',R.SelectedDt, ...
            'loghSpanSigma',P.loghSpanSigma, ...
            'loghYBlockSize',P.loghYBlockSize, ...
            'YGridMode',P.YGridMode,'YGridScale',P.YGridScale, ...
            'minGap',F.minGap,'pinZero',F.pinZero, ...
            'saMaxIter',R.MaxIterSA,'saNStarts',R.NStarts, ...
            'seedInit',R.Seed,'logEvery',0,'historyEvery',P.historyEvery);

        specs{c}=struct('cfg',cfg,'xBaseline',xBaseline(:), ...
            'seed',uint32(R.Seed),'maxEvals',R.MaxEvalsPerStart);
    end
end


function [x,ami]=production_baseline(prod,M,snr,sig)
    if M==32,stage=prod.stages.m32;else,stage=prod.stages.standard;end
    iM=find(stage.axes.M==M,1);
    iS=find(abs(stage.axes.SNR_dB-snr)<1e-12,1);
    iT=find(abs(stage.axes.sigma_X_sq-sig)<1e-12,1);
    if isempty(iM)||isempty(iS)||isempty(iT)
        error('sim_i13d_final_optimizer_benchmark:BaselineLookup','Production baseline not found.');
    end
    z=stage.baseline{iM,iT,iS};
    x=z.x(:);
    if isfield(z,'amiValidated')
        ami=z.amiValidated;
    elseif isfield(z,'AMIValidated')
        ami=z.AMIValidated;
    else
        error('sim_i13d_final_optimizer_benchmark:BaselineAMI','Validated baseline AMI missing.');
    end
end


function sig=build_signature(F,casePlan,sourcePath)
    sig=struct();
    sig.version=F.version;
    sig.sourceProduction=char(sourcePath);
    sig.sourceFreeze=F.sourceProductionFreeze;
    sig.MVec=F.MVec;
    sig.SNRVec=F.SNRVec;
    sig.TurbulenceVec=F.TurbulenceVec;
    sig.Replicates=F.Replicates;
    sig.RestartsByM=F.RestartsByM;
    sig.MaxIterByM=F.MaxIterByM;
    sig.MaxEvalsPerStart=F.MaxEvalsPerStart;
    sig.tieToleranceBits=F.tieToleranceBits;
    sig.minGap=F.minGap;
    sig.pinZero=F.pinZero;
    sig.PatternSearch=F.PatternSearch;
    sig.scheduler=F.scheduler;
    sig.taskGranularity=F.taskGranularity;
    sig.caseIdentities=table2struct(casePlan(:, ...
        {'M','SNRdB','SigmaR2','Replicate','Seed','NStarts','SelectedDt', ...
         'MaxIterSA','MaxEvalsPerStart'}));
end


function [outDir,manifestPath,manifest,resumeRequested]=prepare_run_directory(o,A,F,signatureText)
    resumeRequested=strlength(string(o.ResumeDirectory))>0;

    if resumeRequested
        outDir=char(string(o.ResumeDirectory));
        manifestPath=fullfile(outDir,'i13d_run_manifest.mat');
        if exist(manifestPath,'file')~=2
            error('sim_i13d_final_optimizer_benchmark:ManifestMissing', ...
                'Resume directory has no i13d_run_manifest.mat: %s',outDir);
        end
        X=load(manifestPath,'manifest');
        manifest=X.manifest;
        if ~isfield(manifest,'signatureText') || ~strcmp(manifest.signatureText,signatureText)
            error('sim_i13d_final_optimizer_benchmark:ResumeMismatch', ...
                'Resume directory does not match frozen I.13D signature.');
        end
        fprintf('Resuming I.13D run: %s\n',outDir);
        return;
    end

    if strlength(string(o.OutputDirectory))>0
        parent=char(string(o.OutputDirectory));
    else
        productionRoot=fileparts(A.outputDirectory);
        parent=fullfile(productionRoot,'i13_reviewer12','final_optimizer_benchmark');
    end
    if exist(parent,'dir')~=7,mkdir(parent);end

    token=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    outDir=fullfile(parent,['I13D_' token]);
    [ok,msg]=mkdir(outDir);
    if ~ok,error('sim_i13d_final_optimizer_benchmark:mkdir','%s',msg);end

    manifestPath=fullfile(outDir,'i13d_run_manifest.mat');
    manifest=struct();
    manifest.version=F.version;
    manifest.signatureText=signatureText;
    manifest.createdUTC=utc_now();
    manifest.lastUpdatedUTC=manifest.createdUTC;
    manifest.completed=false;
    manifest.schedulerSegments=repmat(empty_segment(F.ParallelWorkers),0,1);
    manifest.resultPath='';
    save_manifest_atomic(manifestPath,manifest);
end


function [taskResults,segment,trace]=run_dynamic_queue( ...
    taskResults,pending,taskPlan,casePlan,caseSpecs,checkpointPaths, ...
    signatureText,F,o)

    % A fresh process pool is the default because the I.11 production run
    % exposed a stale/long-lived-pool failure mode between stages.
    pp=gcp('nocreate');
    if ~isempty(pp) && o.FreshPool
        delete(pp);
        pp=[];
    end
    if isempty(pp)
        pp=parpool('Processes',double(o.ParallelWorkers));
    elseif pp.NumWorkers~=double(o.ParallelWorkers)
        delete(pp);
        pp=parpool('Processes',double(o.ParallelWorkers));
    end
    if pp.NumWorkers~=double(o.ParallelWorkers)
        error('sim_i13d_final_optimizer_benchmark:PoolSize', ...
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

    % Initial fill.
    while numel(active)<nWorkers && nextPending<=nPending
        rowIdx=pending(nextPending);
        [fut,~]=submit_one(pp,rowIdx,taskPlan,casePlan,caseSpecs,checkpointPaths,signatureText,F);
        active(end+1,1)=fut; %#ok<AGROW>
        activeRows(end+1,1)=rowIdx; %#ok<AGROW>
        nextPending=nextPending+1;
    end

    try
        while ~isempty(active)
            [finishedIdx,res]=fetchNext(active);
            rowIdx=activeRows(finishedIdx);

            active(finishedIdx)=[];
            activeRows(finishedIdx)=[];

            if ~res.success
                cancel(active);
                error('sim_i13d_final_optimizer_benchmark:WorkerFailure', ...
                    'Task %s failed: %s',taskPlan.IdentityKey(rowIdx),res.errorMessage);
            end

            taskResults{rowIdx}=res;
            completed=completed+1;

            % Immediately refill the newly freed worker slot.
            if nextPending<=nPending
                newRow=pending(nextPending);
                [fut,~]=submit_one(pp,newRow,taskPlan,casePlan,caseSpecs,checkpointPaths,signatureText,F);
                active(end+1,1)=fut;
                activeRows(end+1,1)=newRow;
                nextPending=nextPending+1;
            end

            traceElapsed(completed)=toc(wallTimer);
            traceCompleted(completed)=completed;
            traceActive(completed)=numel(active);
            tracePending(completed)=nPending-nextPending+1;
            traceIdentity(completed)=taskPlan.IdentityKey(rowIdx);
            traceWorkerTaskSeconds(completed)=res.taskWallSeconds;

            fprintf('  I13D completed %d/%d this run | total %d/%d | active=%d | pending=%d | %s\n', ...
                completed,nPending,sum(~cellfun(@isempty,taskResults)),height(taskPlan), ...
                numel(active),max(0,nPending-nextPending+1),taskPlan.IdentityKey(rowIdx));
        end
    catch ME
        if ~isempty(active)
            try,cancel(active);catch,end
        end
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
        minActiveBeforeTail=min([nWorkers;trace.ActiveAfterDispatch]);
    end

    segment=struct();
    segment.startedUTC='';
    segment.completedUTC=utc_now();
    segment.workers=nWorkers;
    segment.executedTasks=completed;
    segment.wallSeconds=wallSeconds;
    segment.workerBusySeconds=busySeconds;
    segment.workerCapacitySeconds=capacitySeconds;
    segment.workerUtilization=utilization;
    segment.fullQueueMaintainedUntilTail=fullQueueMaintained;
    segment.minActiveBeforeTail=minActiveBeforeTail;
    segment.scheduler=F.scheduler;

    fprintf('\nI.13D scheduler segment: wall=%.2f h | busy=%.2f worker-h | utilization=%.1f%%\n', ...
        wallSeconds/3600,busySeconds/3600,100*utilization);
    fprintf('Dynamic queue full before tail: %d | min active while pending>0: %d/%d\n', ...
        fullQueueMaintained,minActiveBeforeTail,nWorkers);
end


function [fut,identity]=submit_one(pp,rowIdx,taskPlan,casePlan,caseSpecs,checkpointPaths,signatureText,F)
    t=taskPlan(rowIdx,:);
    c=t.TaskCaseIndex;
    R=casePlan(c,:);
    spec=caseSpecs{c};

    identity=struct();
    identity.IdentityKey=char(t.IdentityKey);
    identity.M=R.M;
    identity.SNRdB=R.SNRdB;
    identity.SigmaR2=R.SigmaR2;
    identity.Replicate=R.Replicate;
    identity.Restart=t.Restart;
    identity.Seed=double(R.Seed);

    fut=parfeval(pp,@i13d_patternsearch_restart_task,1, ...
        spec.cfg,spec.xBaseline,spec.seed,t.Restart,spec.maxEvals, ...
        F.PatternSearch,identity,char(checkpointPaths(rowIdx)),signatureText);
end


function validate_checkpoint(cp,signatureText,identityKey)
    if ~isstruct(cp) || ~isfield(cp,'signatureText') || ~strcmp(cp.signatureText,signatureText)
        error('sim_i13d_final_optimizer_benchmark:CheckpointSignature', ...
            'Checkpoint signature mismatch for %s.',identityKey);
    end
    if ~isfield(cp,'identity') || ~strcmp(string(cp.identity.IdentityKey),string(identityKey))
        error('sim_i13d_final_optimizer_benchmark:CheckpointIdentity', ...
            'Checkpoint identity mismatch for %s.',identityKey);
    end
    if ~isfield(cp,'result')
        error('sim_i13d_final_optimizer_benchmark:CheckpointResult', ...
            'Checkpoint result missing for %s.',identityKey);
    end
end


function [restartSummary,caseSummary,physicalSummary,metrics]=reduce_results(taskResults,taskPlan,casePlan,F)
    n=height(taskPlan);
    M=zeros(n,1); SNRdB=zeros(n,1); SigmaR2=zeros(n,1); Replicate=zeros(n,1);
    Restart=zeros(n,1); Seed=zeros(n,1); InitialFastAMI=zeros(n,1);
    PSFast=zeros(n,1); PSValidated=zeros(n,1); FastValidatorGap=zeros(n,1);
    FunctionEvaluations=zeros(n,1); EvalBudget=zeros(n,1); EvalOverrun=zeros(n,1);
    Iterations=zeros(n,1); Exitflag=zeros(n,1); OptimizationSeconds=zeros(n,1);
    ValidationSeconds=zeros(n,1); TaskWallSeconds=zeros(n,1); IdentityKey=strings(n,1);

    for j=1:n
        r=taskResults{j};
        M(j)=r.M; SNRdB(j)=r.SNRdB; SigmaR2(j)=r.SigmaR2; Replicate(j)=r.Replicate;
        Restart(j)=r.Restart; Seed(j)=r.masterSeed; InitialFastAMI(j)=r.initialFastAMI;
        PSFast(j)=r.bestMIFast; PSValidated(j)=r.bestMIValidated;
        FastValidatorGap(j)=r.validationGap; FunctionEvaluations(j)=r.functionEvaluations;
        EvalBudget(j)=r.maxEvals; EvalOverrun(j)=r.functionEvaluations-r.maxEvals;
        Iterations(j)=r.iterations; Exitflag(j)=r.exitflag;
        OptimizationSeconds(j)=r.optimizationSeconds; ValidationSeconds(j)=r.validationSeconds;
        TaskWallSeconds(j)=r.taskWallSeconds; IdentityKey(j)=r.IdentityKey;
    end

    restartSummary=table(M,SNRdB,SigmaR2,Replicate,Restart,Seed,IdentityKey, ...
        InitialFastAMI,PSFast,PSValidated,FastValidatorGap,FunctionEvaluations, ...
        EvalBudget,EvalOverrun,Iterations,Exitflag,OptimizationSeconds, ...
        ValidationSeconds,TaskWallSeconds);
    restartSummary=sortrows(restartSummary,{'M','SNRdB','SigmaR2','Replicate','Restart'});

    nc=height(casePlan);
    PSValidatedCase=nan(nc,1); PSFastCase=nan(nc,1); PSBestRestart=nan(nc,1);
    PSRestartMean=nan(nc,1); PSRestartStd=nan(nc,1); PSRestartMin=nan(nc,1);
    PSRestartMax=nan(nc,1); PSTotalFunctionEvaluations=nan(nc,1);
    PSMaxEvalOverrun=nan(nc,1); PSTotalTaskSeconds=nan(nc,1);
    PSMaxFastValidatorGap=nan(nc,1); SAMinusPS=nan(nc,1);
    SAGainOverPAM=nan(nc,1); PSGainOverPAM=nan(nc,1);

    for c=1:nc
        C=casePlan(c,:);
        idx=restartSummary.M==C.M & restartSummary.SNRdB==C.SNRdB & ...
            abs(restartSummary.SigmaR2-C.SigmaR2)<1e-12 & ...
            restartSummary.Replicate==C.Replicate;
        R=restartSummary(idx,:);
        if height(R)~=C.NStarts
            error('sim_i13d_final_optimizer_benchmark:ReductionStarts', ...
                'Restart count mismatch for case %d.',c);
        end
        [PSValidatedCase(c),ii]=max(R.PSValidated);
        PSFastCase(c)=R.PSFast(ii);
        PSBestRestart(c)=R.Restart(ii);
        PSRestartMean(c)=mean(R.PSValidated);
        PSRestartStd(c)=std(R.PSValidated,0);
        PSRestartMin(c)=min(R.PSValidated);
        PSRestartMax(c)=max(R.PSValidated);
        PSTotalFunctionEvaluations(c)=sum(R.FunctionEvaluations);
        PSMaxEvalOverrun(c)=max(R.EvalOverrun);
        PSTotalTaskSeconds(c)=sum(R.TaskWallSeconds);
        PSMaxFastValidatorGap(c)=max(abs(R.FastValidatorGap));
        SAMinusPS(c)=C.SAValidated-PSValidatedCase(c);
        SAGainOverPAM(c)=C.SAValidated-C.PAMValidated;
        PSGainOverPAM(c)=PSValidatedCase(c)-C.PAMValidated;
    end

    caseSummary=[casePlan table(PSValidatedCase,PSFastCase,PSBestRestart, ...
        PSRestartMean,PSRestartStd,PSRestartMin,PSRestartMax, ...
        PSTotalFunctionEvaluations,PSMaxEvalOverrun,PSTotalTaskSeconds, ...
        PSMaxFastValidatorGap,SAMinusPS,SAGainOverPAM,PSGainOverPAM)];

    keys=unique(caseSummary(:,{'M','SNRdB','SigmaR2'}),'rows','sorted');
    np=height(keys);
    PAMValidated=nan(np,1); SARep1=nan(np,1); SARep2=nan(np,1);
    PSRep1=nan(np,1); PSRep2=nan(np,1); SAMean=nan(np,1); PSMean=nan(np,1);
    SAStd=nan(np,1); PSStd=nan(np,1); DeltaRep1=nan(np,1); DeltaRep2=nan(np,1);
    SAMinusPSMean=nan(np,1); Outcome=strings(np,1);

    for k=1:np
        idx=caseSummary.M==keys.M(k) & caseSummary.SNRdB==keys.SNRdB(k) & ...
            abs(caseSummary.SigmaR2-keys.SigmaR2(k))<1e-12;
        R=sortrows(caseSummary(idx,:),'Replicate');
        if height(R)~=2 || ~isequal(R.Replicate(:).',[1 2])
            error('sim_i13d_final_optimizer_benchmark:PhysicalReplicates', ...
                'Expected exactly replicates 1 and 2.');
        end
        PAMValidated(k)=mean(R.PAMValidated);
        SARep1(k)=R.SAValidated(1); SARep2(k)=R.SAValidated(2);
        PSRep1(k)=R.PSValidatedCase(1); PSRep2(k)=R.PSValidatedCase(2);
        SAMean(k)=mean([SARep1(k) SARep2(k)]);
        PSMean(k)=mean([PSRep1(k) PSRep2(k)]);
        SAStd(k)=std([SARep1(k) SARep2(k)],0);
        PSStd(k)=std([PSRep1(k) PSRep2(k)],0);
        DeltaRep1(k)=SARep1(k)-PSRep1(k);
        DeltaRep2(k)=SARep2(k)-PSRep2(k);
        SAMinusPSMean(k)=SAMean(k)-PSMean(k);
        if abs(SAMinusPSMean(k))<=F.tieToleranceBits
            Outcome(k)="Equivalent";
        elseif SAMinusPSMean(k)>F.tieToleranceBits
            Outcome(k)="SA better";
        else
            Outcome(k)="Pattern Search better";
        end
    end

    physicalSummary=[keys table(PAMValidated,SARep1,SARep2,PSRep1,PSRep2, ...
        SAMean,PSMean,SAStd,PSStd,DeltaRep1,DeltaRep2,SAMinusPSMean,Outcome)];

    metrics=struct();
    metrics.nRestartTasks=height(restartSummary);
    metrics.nOptimizationCases=height(caseSummary);
    metrics.nPhysicalPoints=height(physicalSummary);
    metrics.tieToleranceBits=F.tieToleranceBits;
    metrics.nEquivalent=sum(physicalSummary.Outcome=="Equivalent");
    metrics.nSABetter=sum(physicalSummary.Outcome=="SA better");
    metrics.nPatternSearchBetter=sum(physicalSummary.Outcome=="Pattern Search better");
    metrics.meanSAMinusPS=mean(physicalSummary.SAMinusPSMean);
    metrics.medianSAMinusPS=median(physicalSummary.SAMinusPSMean);
    metrics.maxSAMinusPS=max(physicalSummary.SAMinusPSMean);
    metrics.minSAMinusPS=min(physicalSummary.SAMinusPSMean);
    metrics.maxAbsDelta=max(abs(physicalSummary.SAMinusPSMean));
    metrics.rmse=sqrt(mean(physicalSummary.SAMinusPSMean.^2));
    metrics.maxPatternSearchFastValidatorGap=max(abs(restartSummary.FastValidatorGap));
    metrics.maxEvaluationOverrun=max(restartSummary.EvalOverrun);
    metrics.totalPatternSearchFunctionEvaluations=sum(restartSummary.FunctionEvaluations);
end


function [P,source]=attach_zhou_if_available(P,outDir)
    source='';
    productionRoot=ancestor_production_root(outDir);
    candidate=fullfile(productionRoot,'i13_reviewer12','zhou_baseline','i13b_zhou_baseline.mat');
    if exist(candidate,'file')~=2
        P.ZhouValidated=nan(height(P),1);
        return;
    end

    Z=load(candidate,'result');
    if ~isfield(Z,'result') || ~isfield(Z.result,'comparison')
        P.ZhouValidated=nan(height(P),1);
        return;
    end

    C=Z.result.comparison;
    z=nan(height(P),1);
    for k=1:height(P)
        idx=C.M==P.M(k) & C.SNRdB==P.SNRdB(k) & abs(C.SigmaX2-P.SigmaR2(k))<1e-12;
        if sum(idx)==1,z(k)=C.ZhouValidated(idx);end
    end
    P.ZhouValidated=z;
    source=candidate;
end


function root=ancestor_production_root(outDir)
    % Expected final benchmark path:
    % production/i13_reviewer12/final_optimizer_benchmark/I13D_...
    root=fileparts(fileparts(fileparts(outDir)));
end


function A=aggregate_scheduler_metrics(segments,defaultWorkers)
    A=struct('nSegments',numel(segments),'wallSeconds',0,'workerBusySeconds',0, ...
        'workerCapacitySeconds',0,'workerUtilization',NaN, ...
        'allSegmentsMaintainedFullQueue',true,'minActiveBeforeTail',defaultWorkers);
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


function s=empty_segment(workers)
    s=struct('startedUTC','','completedUTC','','workers',workers, ...
        'executedTasks',0,'wallSeconds',0,'workerBusySeconds',0, ...
        'workerCapacitySeconds',0,'workerUtilization',NaN, ...
        'fullQueueMaintainedUntilTail',true,'minActiveBeforeTail',workers, ...
        'scheduler','parfeval-dynamic-v1');
end


function print_summary(R)
    M=R.metrics;
    fprintf('\n============================================================\n');
    fprintf('I.13D FINAL BENCHMARK COMPLETE\n');
    fprintf('Physical points: %d | optimization cases: %d | restart tasks: %d\n', ...
        M.nPhysicalPoints,M.nOptimizationCases,M.nRestartTasks);
    fprintf('Equivalence threshold: %.1e bit/symbol\n',M.tieToleranceBits);
    fprintf('Equivalent: %d/%d | SA better: %d/%d | Pattern Search better: %d/%d\n', ...
        M.nEquivalent,M.nPhysicalPoints,M.nSABetter,M.nPhysicalPoints, ...
        M.nPatternSearchBetter,M.nPhysicalPoints);
    fprintf('SA-PS mean/median: %+.6f / %+.6f bit/symbol\n', ...
        M.meanSAMinusPS,M.medianSAMinusPS);
    fprintf('SA-PS min/max: %+.6f / %+.6f | max |delta|=%.6f\n', ...
        M.minSAMinusPS,M.maxSAMinusPS,M.maxAbsDelta);
    fprintf('RMSE: %.6f bit/symbol\n',M.rmse);
    fprintf('Max PS fast-validator gap: %.3e\n',M.maxPatternSearchFastValidatorGap);
    fprintf('Max evaluation-budget overrun: %+g evals\n',M.maxEvaluationOverrun);
    fprintf('PS total function evaluations: %.0f\n',M.totalPatternSearchFunctionEvaluations);

    S=M.scheduler;
    if S.nSegments>0
        fprintf('Scheduler: %s | segments=%d\n',R.config.scheduler,S.nSegments);
        fprintf('Aggregate worker utilization: %.1f%% | busy %.2f worker-h / capacity %.2f worker-h\n', ...
            100*S.workerUtilization,S.workerBusySeconds/3600,S.workerCapacitySeconds/3600);
        fprintf('Full dynamic queue maintained before tail: %d | min active=%d/%d\n', ...
            S.allSegmentsMaintainedFullQueue,S.minActiveBeforeTail,R.config.ParallelWorkers);
    end
    fprintf('Output: %s\n',R.outputDirectory);
    fprintf('============================================================\n');
end


function save_manifest_atomic(path,manifest)
    folder=fileparts(path);
    tmp=[tempname(folder) '.mat'];
    cleanup=onCleanup(@() cleanup_tmp(tmp)); %#ok<NASGU>
    save(tmp,'manifest','-v7');
    [ok,msg]=movefile(tmp,path,'f');
    if ~ok,error('sim_i13d_final_optimizer_benchmark:ManifestSave','%s',msg);end
end

function cleanup_tmp(path)
    if exist(path,'file')==2
        try,delete(path);catch,end
    end
end

function cleanup_pool(pp,doDelete)
    if doDelete && ~isempty(pp)
        try
            if isvalid(pp),delete(pp);end
        catch
        end
    end
end

function s=utc_now()
    s=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0&&mod(v,1)==0;end
