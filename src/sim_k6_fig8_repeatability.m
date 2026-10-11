function R6 = sim_k6_fig8_repeatability(varargin)
%SIM_K6_FIG8_REPEATABILITY  20 independent SA restart trajectories for Fig. 8.
%
% K6 changes only the EXPERIMENT SIZE/SEED relative to frozen I.11.1:
%   M=8, SNR=20 dB, sigma_R^2=0.1; 20 independent Threefry substreams,
%   4000 SA iterations each, same projection, fast AMI and validator.
%
% Source: the exact frozen I.11.1 production bundle. It is read-only.
% The 20 raw restart scores (not multi-start winners) are all used to
% compute sample mean, sample variance (N-1), std and observed range.
%
%   test_k6_fig8_preflight()
%   dry=sim_k6_fig8_repeatability('DryRun',true);
%   R6=sim_k6_fig8_repeatability();
%   R6=sim_k6_fig8_repeatability('ResumeDirectory','/actual/K6_...');
%
% Only new outputs under results/ieee_revision2/commit_k6_fig8.
% No optimizer, evaluator, constraint or production file is modified.
% No global-optimum claim is implied.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'ResumeDirectory','',@text_scalar);
    addParameter(p,'ParallelWorkers',8,@positive_integer);
    addParameter(p,'FreshPool',true,@logical_scalar);
    addParameter(p,'DeletePoolOnFinish',true,@logical_scalar);
    addParameter(p,'MakeFigure',true,@logical_scalar);
    addParameter(p,'DryRun',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    F=i11p1_production_config();
    assert_production_profile(F);

    M=8; snr=20; sigmaR2=0.1; nRestarts=20;
    maxIter=i11p1_iteration_policy(M);
    if maxIter~=4000
        error('K6:Iterations','M=8 must use frozen I.11.1 4000 iterations.');
    end

    prodPath=resolve_production(o.ProductionBundlePath);
    S=load(prodPath,'bundle');
    if ~isfield(S,'bundle') || ~isstruct(S.bundle)
        error('K6:Production','Missing production bundle.');
    end
    B=S.bundle;
    if ~isfield(B,'productionFreeze') || ~isequaln(B.productionFreeze,F)
        error('K6:ProductionFreeze', ...
            'The loaded production scientific freeze differs from I.11.1.');
    end
    if ~isfield(B,'summary') || ~istable(B.summary.caseTable)
        error('K6:ProductionSummary','Missing production case table.');
    end
    C=B.summary.caseTable;
    idx=C.M==M & abs(C.SNRdB-snr)<1e-12 & ...
        abs(C.SigmaX2-sigmaR2)<1e-12;
    sourceCases=sortrows(C(idx,:),'Replicate');
    if height(sourceCases)~=2 || ...
            ~isequal(sourceCases.Replicate(:).',[1 2]) || ...
            any(sourceCases.NStarts~=6) || ...
            any(abs(sourceCases.SelectedDt-F.loghDt)>1e-12)
        error('K6:CaseSource','M=8/20/0.1 original cases missing or changed.');
    end

    stage=B.stages.standard;
    iM=find(stage.axes.M==M);
    iS=find(abs(stage.axes.SNR_dB-snr)<1e-12);
    iT=find(abs(stage.axes.sigma_X_sq-sigmaR2)<1e-12);
    if numel(iM)~=1 || numel(iS)~=1 || numel(iT)~=1
        error('K6:Axes','Invalid original production axes.');
    end
    base=stage.baseline{iM,iT,iS};
    if ~isfield(base,'x') || ~isfield(base,'amiValidated') || ...
            numel(base.x)~=M || ...
            abs(base.amiValidated-mean(sourceCases.AMIValidated-sourceCases.GainBits))>1e-8
        error('K6:Baseline','Production uniform-PAM baseline mismatch.');
    end

    % A NEW fixed case master seed, distinct from both historical production
    % replicate case seeds. One deterministic Threefry substream per restart.
    seed=uint32(fso_result_utils.case_seed(20261011,M,snr,sigmaR2));
    if any(double(seed)==double(sourceCases.Seed))
        error('K6:SeedCollision','New K6 seed duplicates a production case seed.');
    end

    dt=sourceCases.SelectedDt(1);
    xMax=fso_result_utils.feasible_xmax_bound(M,F.P_avg,F.minGap,F.pinZero);
    cfg=build_fso_config(M,F.P_avg,snr,sigmaR2, ...
        'FastFadingMethod',F.FastFadingMethod, ...
        'loghDt',dt, ...
        'loghSpanSigma',F.loghSpanSigma, ...
        'loghYBlockSize',F.loghYBlockSize, ...
        'YGridMode',F.YGridMode, 'YGridScale',F.YGridScale, ...
        'xMaxBound',xMax, ...
        'saMaxIter',maxIter,'saNStarts',nRestarts, ...
        'saUseParallel',false, ...
        'minGap',F.minGap,'pinZero',F.pinZero, ...
        'seedInit',double(seed),'logEvery',0, ...
        'historyEvery',F.historyEvery);
    assert_exact_config(cfg,F,M,snr,sigmaR2,dt,maxIter,nRestarts,seed);
    AMI_functions.assert_constellation_feasible(base.x(:),cfg,1e-10);

    Restart=(1:nRestarts).';
    MasterSeed=repmat(double(seed),nRestarts,1);
    Substream=Restart;
    IdentityKey=arrayfun(@(r)string(sprintf( ...
        'K6_M8_SNR20_sig100_r%02d',r)),Restart);
    taskPlan=table(Restart,MasterSeed,Substream,IdentityKey);
    if numel(unique(taskPlan.IdentityKey))~=nRestarts
        error('K6:TaskIdentity','Duplicate restart identity.');
    end
    signature=struct('version','K6-FIG8-I11p1-v1', ...
        'sourceProduction',canonical(prodPath), ...
        'productionFreeze',F, ...
        'M',M,'SNRdB',snr,'SigmaR2',sigmaR2, ...
        'nRestarts',nRestarts,'maxIter',maxIter, ...
        'seedPolicy','new Threefry master seed, unique substreams 1:20', ...
        'masterSeed',double(seed), ...
        'cfg',fso_result_utils.config_snapshot(cfg), ...
        'taskIds',{cellstr(taskPlan.IdentityKey)});
    signatureText=char(jsonencode(signature));

    fprintf('\n============================================================\n');
    fprintf('COMMIT K6 - FIG. 8 SA REPEATABILITY (%s)\n',signature.version);
    fprintf('M=%d | SNR=%g dB | sigma_R^2=%.1f | P_avg=%.1f\n', ...
        M,snr,sigmaR2,F.P_avg);
    fprintf('d_min=%.3f | first level pinned to zero | no xMax shaping bound\n',F.minGap);
    fprintf('SA: %d independent restarts x %d iterations = %d iterations\n', ...
        nRestarts,maxIter,nRestarts*maxIter);
    fprintf('T0=%.4f | Tf=%.4g | BaseStd0=%.3f | block=%d\n', ...
        F.T0,F.Tf,F.BaseStd0,F.ItersPerTemp);
    fprintf('Evaluator=logh | candidate-adaptive y-grid | logh dt=%.3g\n',dt);
    fprintf('Threefry masterSeed=%u | substreams=1:20\n',seed);
    fprintf('Source production (read only): %s\n',prodPath);
    fprintf('No global-optimum claim is implied.\n');
    fprintf('============================================================\n');

    if o.DryRun
        R6=struct('dryRun',true,'version',signature.version, ...
            'config',fso_result_utils.config_snapshot(cfg), ...
            'sourceProduction',prodPath,'sourceCases',sourceCases, ...
            'taskPlan',taskPlan,'signature',signature, ...
            'nRestarts',nRestarts,'totalSAIterations',nRestarts*maxIter, ...
            'outputDirectory','');
        return;
    end

    outDir=prepare_run_directory(o);
    cpDir=fullfile(outDir,'checkpoints');
    if exist(cpDir,'dir')~=7,mkdir(cpDir);end
    manifestPath=fullfile(outDir,'K6_run_manifest.mat');
    if exist(manifestPath,'file')==2
        Z=load(manifestPath,'signatureText');
        if ~isfield(Z,'signatureText') || ~strcmp(Z.signatureText,signatureText)
            error('K6:ResumeMismatch','Saved K6 run has a different scientific signature.');
        end
    else
        if strlength(string(o.ResumeDirectory))>0
            error('K6:ResumeManifest','Resume directory has no K6 manifest.');
        end
        save(manifestPath,'signature','signatureText');
    end
    fprintf('K6 run directory: %s\n',outDir);

    paths=cell(nRestarts,1);
    identities=cell(nRestarts,1);
    restartResults=cell(nRestarts,1);
    reused=false(nRestarts,1);

    for k=1:nRestarts
        identities{k}=make_identity(M,snr,sigmaR2,seed,k,taskPlan.IdentityKey(k));
        paths{k}=fullfile(cpDir,[char(taskPlan.IdentityKey(k)) '.mat']);
        if exist(paths{k},'file')~=2,continue;end
        z=load(paths{k},'restartResult','checkpointMeta');
        if ~isfield(z,'checkpointMeta') || ~isfield(z,'restartResult') || ...
                ~strcmp(z.checkpointMeta.signatureText,signatureText) || ...
                ~isequaln(z.checkpointMeta.identity,identities{k})
            error('K6:CheckpointMismatch','Invalid prior checkpoint: %s',paths{k});
        end
        validate_restart(z.restartResult,cfg,seed,k);
        restartResults{k}=z.restartResult;
        reused(k)=true;
    end
    pending=find(~reused);
    fprintf('Validated checkpoint reuse: %d/%d | pending %d\n', ...
        nnz(reused),nRestarts,numel(pending));

    runStats=struct('workers',0,'wallSeconds',0,'workerBusySeconds',0, ...
        'workerUtilization',NaN,'completedTasks',0, ...
        'checkpointsReused',nnz(reused));
    if ~isempty(pending)
        if o.ParallelWorkers==1
            wall=tic;
            busy=0;
            for i=1:numel(pending)
                k=pending(i);
                r=sa_restart_task_checkpointed(cfg,base.x,seed,k, ...
                    sprintf('[K6|%02d]',k),paths{k},signatureText,identities{k});
                validate_restart(r,cfg,seed,k);
                restartResults{k}=r;
                busy=busy+r.taskWallSeconds;
                fprintf('K6 DONE %d/%d | restart=%02d | validated=%.8f\n', ...
                    i,numel(pending),k,r.bestMIValidated);
            end
            runStats.workers=1;
            runStats.wallSeconds=toc(wall);
            runStats.workerBusySeconds=busy;
            runStats.workerUtilization=busy/max(runStats.wallSeconds,eps);
        else
            pp=gcp('nocreate');
            if ~isempty(pp) && o.FreshPool,delete(pp);pp=[];end
            if isempty(pp)
                pp=parpool('Processes',o.ParallelWorkers);
            elseif pp.NumWorkers~=o.ParallelWorkers
                delete(pp);pp=parpool('Processes',o.ParallelWorkers);
            end
            if o.DeletePoolOnFinish
                poolCleanup=onCleanup(@() delete_pool(pp)); %#ok<NASGU>
            end
            runStats.workers=pp.NumWorkers;
            active=parallel.FevalFuture.empty(0,1);
            activeRestart=zeros(0,1);
            next=1;done=0;busy=0;
            wall=tic;
            try
                while numel(active)<pp.NumWorkers && next<=numel(pending)
                    k=pending(next);
                    active(end+1,1)=submit_one(pp,cfg,base.x,seed, ...
                        k,paths{k},signatureText,identities{k}); %#ok<AGROW>
                    activeRestart(end+1,1)=k; %#ok<AGROW>
                    next=next+1;
                end
                while ~isempty(active)
                    [pos,r]=fetchNext(active);
                    k=activeRestart(pos);
                    active(pos)=[];activeRestart(pos)=[];
                    validate_restart(r,cfg,seed,k);
                    restartResults{k}=r;done=done+1;
                    busy=busy+r.taskWallSeconds;
                    if next<=numel(pending)
                        kNext=pending(next);
                        active(end+1,1)=submit_one(pp,cfg,base.x,seed, ...
                            kNext,paths{kNext},signatureText, ...
                            identities{kNext}); %#ok<AGROW>
                        activeRestart(end+1,1)=kNext; %#ok<AGROW>
                        next=next+1;
                    end
                    fprintf(['K6 DONE %2d/%2d | restart=%02d | ' ...
                        'validated=%.8f | active=%d | unsubmitted=%d\n'], ...
                        done,numel(pending),k,r.bestMIValidated, ...
                        numel(active),numel(pending)-next+1);
                end
            catch ME
                if ~isempty(active),try,cancel(active);catch,end,end
                rethrow(ME);
            end
            runStats.wallSeconds=toc(wall);
            runStats.workerBusySeconds=busy;
            runStats.workerUtilization=busy/max(pp.NumWorkers*runStats.wallSeconds,eps);
        end
        runStats.completedTasks=numel(pending);
    end

    % Assemble in RESTART INDEX ORDER. No selection/censoring of trajectories.
    restartTable=make_restart_table(restartResults,cfg,seed,nRestarts);
    constellationTable=make_constellation_table(restartResults,M);
    [statistics,summaryTable]=summarize_all(restartTable);
    if height(restartTable)~=nRestarts || height(constellationTable)~=nRestarts
        error('K6:ResultCount','Expected exactly 20 complete restart rows.');
    end
    if any(~isfinite(restartTable.ValidatedAMI))
        error('K6:ValidatedAMI','Every SA trajectory must be validated.');
    end

    writetable(restartTable,fullfile(outDir,'K6_restart_results.csv'));
    writetable(constellationTable,fullfile(outDir,'K6_constellation_levels.csv'));
    writetable(summaryTable,fullfile(outDir,'K6_repeatability_summary.csv'));

    figurePaths=struct();
    if o.MakeFigure
        figurePaths=plot_k6_fig8_repeatability(restartTable,statistics, ...
            'OutputDirectory',outDir);
    end

    R6=struct();
    R6.version=signature.version;
    R6.sourceProduction=prodPath;
    R6.productionFreeze=F;
    R6.config=fso_result_utils.config_snapshot(cfg);
    R6.signature=signature;
    R6.taskPlan=taskPlan;
    R6.restartResults=restartResults;
    R6.restartTable=restartTable;
    R6.constellationTable=constellationTable;
    R6.statistics=statistics;
    R6.summaryTable=summaryTable;
    R6.runStats=runStats;
    R6.figurePaths=figurePaths;
    R6.outputDirectory=outDir;
    R6.globalOptimumClaim=false;
    save(fullfile(outDir,'K6_repeatability_results.mat'),'R6','-v7.3');

    fprintf('\n============================================================\n');
    fprintf('COMMIT K6 COMPLETE - %d independently validated SA trajectories\n',nRestarts);
    fprintf('Mean %.9f | sample variance %.9g | sample std %.9g bits/symbol\n', ...
        statistics.MeanAMI,statistics.SampleVarianceAMI,statistics.SampleStdAMI);
    fprintf('Min %.9f | max %.9f | range %.9g | median %.9f\n', ...
        statistics.MinAMI,statistics.MaxAMI,statistics.RangeAMI,statistics.MedianAMI);
    fprintf('Max |fast-validator| = %.3e bit/symbol\n', ...
        statistics.MaxAbsFastValidatorGap);
    fprintf('Wall %.2f min | worker utilization %.1f%% (executed tasks)\n', ...
        runStats.wallSeconds/60,100*runStats.workerUtilization);
    fprintf('Output: %s\n',outDir);
    fprintf('No multi-start winner averaging; all 20 raw restart results included.\n');
    fprintf('No global-optimum claim is made.\n');
    fprintf('============================================================\n');
    disp(summaryTable);
end


function assert_production_profile(F)
    if ~strcmp(F.freezeVersion,'I11.1-v1') || ...
            ~isequal(F.MVec,[4 8 16 32]) || ...
            ~isequal(F.SNRVec,5:5:30) || ...
            ~isequal(F.TurbulenceVec,[0 0.1 0.2 0.3]) || ...
            F.P_avg~=1 || F.minGap~=0.01 || ~F.pinZero || ...
            ~isequal(F.RestartsByM,[6 6 8 8]) || ...
            ~isequal(F.MaxIterByM,[4000 4000 4000 8000]) || ...
            F.nReplicates~=2 || F.T0~=0.4 || F.Tf~=1e-3 || ...
            F.BaseStd0~=0.15 || F.ItersPerTemp~=50 || ...
            ~strcmp(F.FastFadingMethod,'logh') || ...
            ~strcmp(F.LoghDtPolicy,'case-adaptive') || ...
            F.loghDt~=0.01 || F.loghDtRefined~=0.005 || ...
            F.loghSpanSigma~=8 || F.loghYBlockSize~=512 || ...
            ~strcmp(F.YGridMode,'candidate-adaptive') || ...
            F.YGridScale~=1.1 || F.historyEvery~=50 || ...
            F.baseSeed~=20260928
        error('K6:FrozenConfig', ...
            'Frozen I.11.1 scientific policy changed; aborting K6.');
    end
end


function assert_exact_config(cfg,F,M,snr,sig,dt,maxIter,nRestarts,seed)
    sa=cfg.SA;
    if cfg.M~=M || cfg.P_avg~=F.P_avg || cfg.SNR_dB~=snr || ...
            abs(cfg.sigma_X_sq-sig)>1e-12 || cfg.R~=1 || ...
            abs(cfg.sigma_n_sq-1e-2)>1e-12 || ...
            ~isequal(cfg.px,ones(M,1)/M) || ...
            ~strcmp(cfg.fastFadingMethod,F.FastFadingMethod) || ...
            ~strcmp(cfg.yGridMode,F.YGridMode) || ...
            cfg.yGridScale~=F.YGridScale || ...
            cfg.loghDt~=dt || cfg.loghSpanSigma~=F.loghSpanSigma || ...
            cfg.loghYBlockSize~=F.loghYBlockSize || ...
            cfg.xMaxBound~=fso_result_utils.feasible_xmax_bound( ...
                M,F.P_avg,F.minGap,F.pinZero) || ...
            sa.maxIter~=maxIter || sa.nStarts~=nRestarts || ...
            sa.T0~=F.T0 || sa.Tf~=F.Tf || ...
            sa.baseStd0~=F.BaseStd0 || sa.itersPerTemp~=F.ItersPerTemp || ...
            abs(sa.coolingRate-exp(log(F.Tf/F.T0)/ ...
                ceil(maxIter/F.ItersPerTemp)))>1e-12 || ...
            sa.minGap~=F.minGap || ~sa.pinZero || ...
            ~sa.enforce_sort || ~sa.enforce_power || ...
            sa.powerConstraint~="mean" || ...
            ~sa.imdd_mode || sa.seedInit~=double(seed) || ...
            sa.historyEvery~=F.historyEvery || sa.logEvery~=0 || ...
            sa.useParallel || sa.baseStdMin~=1e-3 || ...
            sa.baseStdMax~=0.5 || sa.targetAccLo~=0.2 || ...
            sa.targetAccHi~=0.6 || sa.baseStdGrow~=1.25 || ...
            sa.baseStdShrink~=0.8
        error('K6:ConfigMismatch', ...
            'SA, objective, power, spacing, or initialization differs from I.11.1.');
    end
end


function f=submit_one(pp,cfg,baseline,seed,k,path,sig,identity)
    f=parfeval(pp,@sa_restart_task_checkpointed,1, ...
        cfg,baseline,seed,k,sprintf('[K6|%02d]',k),path,sig,identity);
end


function validate_restart(r,cfg,seed,k)
    required={'startIndex','masterSeed','substream','bestMIValidated', ...
        'bestMIFast','validationGap','x_best','nIterations', ...
        'runtime','taskWallSeconds','validationRuntime', ...
        'nAccepted','acceptanceRate','initialMI','history'};
    if ~isstruct(r)||~all(isfield(r,required))
        error('K6:IncompleteRestart','Malformed SA restart %d.',k);
    end
    if r.startIndex~=k || double(r.masterSeed)~=double(seed) || ...
            r.substream~=k || r.nIterations~=4000 || ...
            ~isfinite(r.bestMIValidated) || ~isfinite(r.bestMIFast) || ...
            ~isfinite(r.initialMI) || ~isfinite(r.validationGap) || ...
            abs(r.validationGap-(r.bestMIFast-r.bestMIValidated))>1e-10 || ...
            ~isfinite(r.runtime) || ~isfinite(r.validationRuntime) || ...
            ~isfinite(r.taskWallSeconds) || ...
            r.acceptanceRate<0 || r.acceptanceRate>1
        error('K6:RestartValue','Invalid validated SA restart %d.',k);
    end
    AMI_functions.assert_constellation_feasible(r.x_best,cfg,1e-10);
    if abs(mean(r.x_best)-cfg.P_avg)>1e-9 || ...
            numel(r.x_best)~=8 || ...
            min(diff(r.x_best(:)))<cfg.SA.minGap-1e-9
        error('K6:Geometry','Restart %d constellation violates constraints.',k);
    end
    if numel(r.history.iter)<2 || r.history.iter(1)~=0 || ...
            r.history.iter(end)~=4000
        error('K6:History','Incomplete 4000-iteration convergence history.');
    end
end


function T=make_restart_table(results,cfg,seed,n)
    Restart=(1:n).';
    MasterSeed=repmat(double(seed),n,1);
    Substream=Restart;
    InitialFastAMI=zeros(n,1);BestFastAMI=zeros(n,1);
    ValidatedAMI=zeros(n,1);FastValidatorGap=zeros(n,1);
    SASeconds=zeros(n,1);ValidationSeconds=zeros(n,1);
    TaskWallSeconds=zeros(n,1);NIterations=zeros(n,1);
    AcceptedMoves=zeros(n,1);WorseAcceptedMoves=zeros(n,1);
    AcceptanceRate=zeros(n,1);MinAdjacentGap=zeros(n,1);
    MeanIntensity=zeros(n,1);MaxIntensity=zeros(n,1);
    K95=zeros(n,1);K99=zeros(n,1);
    for k=1:n
        r=results{k};
        validate_restart(r,cfg,seed,k);
        InitialFastAMI(k)=r.initialMI;
        BestFastAMI(k)=r.bestMIFast;
        ValidatedAMI(k)=r.bestMIValidated;
        FastValidatorGap(k)=r.validationGap;
        SASeconds(k)=r.runtime;
        ValidationSeconds(k)=r.validationRuntime;
        TaskWallSeconds(k)=r.taskWallSeconds;
        NIterations(k)=r.nIterations;
        AcceptedMoves(k)=r.nAccepted;
        WorseAcceptedMoves(k)=r.nWorseAccepted;
        AcceptanceRate(k)=r.acceptanceRate;
        MinAdjacentGap(k)=min(diff(r.x_best(:)));
        MeanIntensity(k)=mean(r.x_best);
        MaxIntensity(k)=max(r.x_best);
        K95(k)=convergence_iter(r.history,r.initialMI,r.bestMIFast,.95);
        K99(k)=convergence_iter(r.history,r.initialMI,r.bestMIFast,.99);
    end
    T=table(Restart,MasterSeed,Substream,InitialFastAMI, ...
        BestFastAMI,ValidatedAMI,FastValidatorGap,NIterations, ...
        AcceptedMoves,WorseAcceptedMoves,AcceptanceRate, ...
        SASeconds,ValidationSeconds,TaskWallSeconds, ...
        MinAdjacentGap,MeanIntensity,MaxIntensity,K95,K99);
end


function k=convergence_iter(h,initial,best,fraction)
    improvement=best-initial;
    if ~isfinite(improvement)||improvement<=0
        k=0;return;
    end
    pos=find(h.bestMI>=initial+fraction*improvement,1,'first');
    if isempty(pos),k=h.iter(end);else,k=h.iter(pos);end
end


function T=make_constellation_table(results,M)
    n=numel(results);
    levels=zeros(n,M);ValidatedAMI=zeros(n,1);
    for k=1:n
        x=results{k}.x_best(:).';
        levels(k,:)=x;
        ValidatedAMI(k)=results{k}.bestMIValidated;
    end
    T=table((1:n).',ValidatedAMI,'VariableNames',{'Restart','ValidatedAMI'});
    for j=1:M
        T.(sprintf('x%d',j))=levels(:,j);
    end
end


function [stats,S]=summarize_all(T)
    vals=T.ValidatedAMI;
    if height(T)~=20 || any(~isfinite(vals))
        error('K6:Stats','Expected exactly 20 finite validated results.');
    end
    [low,worst]=min(vals);[high,best]=max(vals);
    stats=struct();
    stats.N=height(T);
    stats.MeanAMI=mean(vals);
    stats.SampleVarianceAMI=var(vals,0);
    stats.SampleStdAMI=std(vals,0);
    stats.MedianAMI=median(vals);
    stats.MinAMI=low;
    stats.MaxAMI=high;
    stats.RangeAMI=high-low;
    stats.BestRestart=T.Restart(best);
    stats.WorstRestart=T.Restart(worst);
    stats.MaxAbsFastValidatorGap=max(abs(T.FastValidatorGap));
    stats.MeanAbsFastValidatorGap=mean(abs(T.FastValidatorGap));
    stats.MeanSASeconds=mean(T.SASeconds);
    stats.TotalSASeconds=sum(T.SASeconds);
    stats.MeanAcceptanceRate=mean(T.AcceptanceRate);
    S=struct2table(stats);
end


function path=resolve_production(given)
    if strlength(string(given))>0
        path=char(string(given));
    else
        root=fullfile(fso_result_utils.default_results_root(),'i11p1_production');
        D=dir(fullfile(root,'*','i11p1ProductionBundle.mat'));
        if numel(D)~=1
            error('K6:ProductionDiscovery', ...
                ['Found %d production bundles. ' ...
                 'Supply ProductionBundlePath explicitly.'],numel(D));
        end
        path=fullfile(D(1).folder,D(1).name);
    end
    if exist(path,'file')~=2
        error('K6:ProductionMissing','Missing production MAT: %s',path);
    end
end


function outDir=prepare_run_directory(o)
    root=canonical(fso_result_utils.revision2_results_root());
    if strlength(string(o.ResumeDirectory))>0
        outDir=char(string(o.ResumeDirectory));
        if exist(outDir,'dir')~=7
            error('K6:ResumeDirectory','Resume directory does not exist.');
        end
    elseif strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
        if exist(outDir,'dir')==7
            error('K6:OutputExists','Use ResumeDirectory for an existing K6 folder.');
        end
    else
        parent=fullfile(root,'commit_k6_fig8');
        if exist(parent,'dir')~=7,mkdir(parent);end
        token=char(datetime('now','TimeZone','UTC', ...
            'Format','yyyyMMdd''T''HHmmssSSS''Z'''));
        outDir=fullfile(parent,['K6_' token]);
    end
    target=canonical(outDir);
    if ~(strcmp(target,root)||startsWith(target,[root filesep]))
        error('K6:ResultRoot','Output must be under results/ieee_revision2.');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('K6:mkdir','%s',msg);end
    end
end


function id=make_identity(M,snr,sig,seed,k,key)
    id=struct('IdentityKey',char(key),'M',M,'SNRdB',snr, ...
        'SigmaR2',sig,'MasterSeed',double(seed),'Restart',k);
end


function delete_pool(pp)
    if ~isempty(pp)
        try,delete(pp);catch,end
    end
end
function s=canonical(path)
    s=char(java.io.File(char(string(path))).getCanonicalPath());
end
function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_integer(v)
    tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1&&v==floor(v);
end
