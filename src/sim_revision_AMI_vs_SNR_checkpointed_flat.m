function bundle = sim_revision_AMI_vs_SNR_checkpointed_flat(varargin)
%SIM_REVISION_AMI_VS_SNR_CHECKPOINTED_FLAT  I.11-ready flat restart scheduler.
%
% The atomic work item is one independently validated SA restart:
%   (M,SNR_dB,sigma_X^2,replicate,restart).
% All tasks share one process-pool queue; there is no nested parallelism.
%
% I.11 additions over I.10.5A:
%   - restart-count policy may depend on M (I11: 6 for M<=8, 8 for M>=16);
%   - every completed restart can be checkpointed atomically;
%   - ResumeDirectory reuses only checkpoints whose scientific run signature
%     and task identity match exactly;
%   - Uniform-PAM baselines are cached as well;
%   - worker count and ordering are operational controls and may change on
%     resume without changing the stochastic trajectories.
%
% No global-optimum claim is implied by the restart policy.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'MVec',[4 8 16 32],@valid_mvec);
    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'SNRVec',5:5:30,@numeric_vector);
    addParameter(p,'TurbulenceVec',[0 0.1 0.2 0.3],@nonnegative_vector);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'pinZero',true,@logical_scalar);
    addParameter(p,'FastFadingMethod','logh',@text_scalar);
    addParameter(p,'ghN_h',400,@positive_integer);
    addParameter(p,'LoghDtPolicy','case-adaptive',@text_scalar);
    addParameter(p,'loghDt',0.01,@positive_scalar);
    addParameter(p,'loghDtRefined',0.005,@positive_scalar);
    addParameter(p,'loghSpanSigma',8,@positive_scalar);
    addParameter(p,'loghYBlockSize',512,@positive_integer);
    addParameter(p,'YGridMode','candidate-adaptive',@text_scalar);
    addParameter(p,'YGridScale',1.1,@(v) positive_scalar(v)&&v>=1);

    addParameter(p,'saMaxIter',4000,@positive_integer);
    addParameter(p,'RestartPolicy','fixed',@text_scalar);
    addParameter(p,'saNStarts',6,@positive_integer);
    addParameter(p,'nReplicates',2,@positive_integer);
    addParameter(p,'T0',0.4,@positive_scalar);
    addParameter(p,'Tf',1e-3,@positive_scalar);
    addParameter(p,'BaseStd0',0.15,@positive_scalar);
    addParameter(p,'ItersPerTemp',50,@positive_integer);
    addParameter(p,'baseSeed',20260928,@nonnegative_scalar);
    addParameter(p,'historyEvery',50,@positive_integer);

    addParameter(p,'UseParallel',true,@logical_scalar);
    addParameter(p,'ParallelWorkers',[],@(v) isempty(v)||positive_integer(v));
    addParameter(p,'LongestFirst',true,@logical_scalar);
    addParameter(p,'Checkpoint',true,@logical_scalar);
    addParameter(p,'ResumeDirectory','',@text_scalar);
    addParameter(p,'WriteCSV',true,@logical_scalar);
    addParameter(p,'RunLabel','',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    parse(p,varargin{:}); o=p.Results;

    if ~(o.T0>o.Tf)
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:TemperatureOrder','T0 must exceed Tf.');
    end
    if o.loghDtRefined>o.loghDt
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:DtOrder','Refined dt must not exceed base dt.');
    end

    fastMethod=lower(string(o.FastFadingMethod));
    yMode=lower(string(o.YGridMode));
    dtPolicy=lower(string(o.LoghDtPolicy));
    restartPolicy=lower(string(o.RestartPolicy));
    if ~ismember(fastMethod,["logh","gh"])
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:FastMethod','FastFadingMethod must be logh or gh.');
    end
    if ~ismember(yMode,["fixed-bound","candidate-adaptive"])
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:YGridMode','Unsupported YGridMode.');
    end
    if ~ismember(dtPolicy,["case-adaptive","fixed"])
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:DtPolicy','Unsupported LoghDtPolicy.');
    end
    if ~ismember(restartPolicy,["fixed","i11"])
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:RestartPolicy', ...
            'RestartPolicy must be ''fixed'' or ''i11''.');
    end
    if fastMethod=="gh" && yMode=="candidate-adaptive"
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:AdaptiveRequiresLogh', ...
            'candidate-adaptive y-grid is supported only by the log-h evaluator.');
    end
    if strlength(string(o.ResumeDirectory))>0 && ~o.Checkpoint
        error('sim_revision_AMI_vs_SNR_checkpointed_flat:ResumeRequiresCheckpoint', ...
            'ResumeDirectory requires Checkpoint=true.');
    end
    o.FastFadingMethod=char(fastMethod); o.YGridMode=char(yMode);
    o.LoghDtPolicy=char(dtPolicy); o.RestartPolicy=char(restartPolicy);

    MVec=unique(double(o.MVec(:).'),'stable');
    snrVec=unique(double(o.SNRVec(:).'),'stable');
    sigVec=unique(double(o.TurbulenceVec(:).'),'stable');
    Pavg=double(o.P_avg); nRep=double(o.nReplicates);
    validate_spacing_feasibility(MVec,Pavg,double(o.minGap));

    startsByM=zeros(size(MVec)); restartMeta=cell(size(MVec));
    for k=1:numel(MVec)
        if restartPolicy=="i11"
            [startsByM(k),restartMeta{k}]=i11_restart_policy(MVec(k));
        else
            startsByM(k)=double(o.saNStarts);
            restartMeta{k}=struct('ruleVersion','fixed','M',MVec(k), ...
                'nStarts',startsByM(k),'reason','Explicit fixed restart count');
        end
    end

    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec);
    nPhysical=nM*nSig*nSNR; nCases=nPhysical*nRep;
    nTasks=numel(sigVec)*numel(snrVec)*nRep*sum(startsByM);

    xMaxBounds=nan(1,nM);
    for iM=1:nM
        xMaxBounds(iM)=fso_result_utils.feasible_xmax_bound( ...
            MVec(iM),Pavg,double(o.minGap),logical(o.pinZero));
    end
    [dtGrid,dtMeta]=resolve_dt_grid(MVec,sigVec,snrVec,o,fastMethod);
    refinedMask=cellfun(@(m) isfield(m,'isRefined')&&m.isRefined,dtMeta);

    signature=build_signature(MVec,snrVec,sigVec,Pavg,startsByM,o);
    signatureText=char(jsonencode(signature));

    resumeRequested=strlength(string(o.ResumeDirectory))>0;
    if resumeRequested
        outDir=char(string(o.ResumeDirectory));
        if exist(outDir,'dir')~=7
            error('sim_revision_AMI_vs_SNR_checkpointed_flat:ResumeDirectoryMissing', ...
                'ResumeDirectory does not exist: %s',outDir);
        end
        [~,runId]=fileparts(outDir);
    else
        parentDir=fullfile(char(o.ResultsRoot),'revised_ami_vs_snr_checkpointed_flat');
        label=sanitize_label(o.RunLabel);
        token=sprintf('I%d_P%s_R%d',o.saMaxIter,char(restartPolicy),o.nReplicates);
        utc=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
        [outDir,runId]=make_unique_dir(parentDir,sprintf('%s_%s_%s',label,token,utc));
    end

    caseDir=fullfile(outDir,'cases'); if ~exist(caseDir,'dir'),mkdir(caseDir);end
    checkpointDir=fullfile(outDir,'checkpoints'); if o.Checkpoint && ~exist(checkpointDir,'dir'),mkdir(checkpointDir);end
    manifestPath=fullfile(outDir,'run_manifest.mat');

    if resumeRequested
        if exist(manifestPath,'file')~=2
            error('sim_revision_AMI_vs_SNR_checkpointed_flat:ResumeManifestMissing', ...
                'Resume directory has no run_manifest.mat.');
        end
        S=load(manifestPath,'manifest');
        if ~isfield(S,'manifest') || ~isfield(S.manifest,'signatureText') || ...
                ~strcmp(S.manifest.signatureText,signatureText)
            error('sim_revision_AMI_vs_SNR_checkpointed_flat:ResumeConfigMismatch', ...
                ['Resume configuration does not match the original scientific signature. ' ...
                 'Do not reuse checkpoints across different scientific settings.']);
        end
        manifest=S.manifest;
    else
        manifest=struct();
        manifest.schedulerVersion='I11-checkpointed-flat-v1';
        manifest.createdUTC=char(datetime('now','TimeZone','UTC', ...
            'Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
        manifest.signature=signature;
        manifest.signatureText=signatureText;
        manifest.runId=runId;
        manifest.runLabel=char(string(o.RunLabel));
        save(manifestPath,'manifest','-v7');
    end

    runMeta=fso_result_utils.run_metadata(mfilename);
    runMeta.runId=runId; runMeta.runLabel=char(string(o.RunLabel));
    runMeta.outputDirectory=outDir; runMeta.resumeRequested=resumeRequested;

    saProfile=struct('T0',double(o.T0),'Tf',double(o.Tf), ...
        'baseStd0',double(o.BaseStd0),'itersPerTemp',double(o.ItersPerTemp));

    fprintf('\n============================================================\n');
    fprintf('I.11 checkpointed flat restart-task scheduler\n');
    fprintf('Run ID: %s | resume=%d\n',runId,resumeRequested);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | replicates=%d\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),nRep);
    fprintf('Restart policy=%s | starts by M=%s | maxIter=%d\n', ...
        char(restartPolicy),mat2str(startsByM),o.saMaxIter);
    fprintf('Flattened restart tasks=%d | longest-first=%d | checkpoint=%d\n', ...
        nTasks,logical(o.LongestFirst),logical(o.Checkpoint));
    fprintf('Refined physical points: %d/%d\n',sum(refinedMask(:)),nPhysical);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    baselinePath=fullfile(outDir,'baseline_cache.mat'); baselineReused=false;
    if o.Checkpoint && exist(baselinePath,'file')==2
        try
            B=load(baselinePath,'baseline','baselineSignatureText');
            if isfield(B,'baseline') && isfield(B,'baselineSignatureText') && ...
                    strcmp(B.baselineSignatureText,signatureText)
                baseline=B.baseline; baselineReused=true;
            else
                baseline=[];
            end
        catch
            baseline=[];
        end
    else
        baseline=[];
    end
    if isempty(baseline)
        baseline=build_baselines(MVec,sigVec,snrVec,Pavg,o,xMaxBounds,dtGrid,dtMeta,saProfile,startsByM);
        if o.Checkpoint
            baselineSignatureText=signatureText; %#ok<NASGU>
            tmp=[tempname(outDir) '.mat']; cleanup=onCleanup(@() cleanup_tmp(tmp)); %#ok<NASGU>
            save(tmp,'baseline','baselineSignatureText','-v7');
            [ok,msg]=movefile(tmp,baselinePath,'f');
            if ~ok,error('sim_revision_AMI_vs_SNR_checkpointed_flat:BaselineSave','%s',msg);end
        end
    else
        fprintf('Reusing cached Uniform-PAM baselines.\n');
    end

    caseSpecs=cell(nCases,1); q=0;
    for iM=1:nM
        for iSig=1:nSig
            for iSNR=1:nSNR
                for iR=1:nRep
                    q=q+1;
                    repBaseSeed=double(o.baseSeed)+(iR-1)*10000019;
                    seed=fso_result_utils.case_seed(repBaseSeed,MVec(iM),snrVec(iSNR),sigVec(iSig));
                    s=struct();
                    s.caseId=q; s.iM=iM; s.iSig=iSig; s.iSNR=iSNR; s.iR=iR;
                    s.M=MVec(iM); s.SNRdB=snrVec(iSNR); s.sigmaX2=sigVec(iSig);
                    s.seed=seed; s.selectedDt=dtGrid(iM,iSig,iSNR); s.dtMeta=dtMeta{iM,iSig,iSNR};
                    s.xMaxBound=xMaxBounds(iM); s.baseline=baseline{iM,iSig,iSNR};
                    s.nStarts=startsByM(iM); s.restartMeta=restartMeta{iM};
                    s.label=sprintf('[I11|M%d|SNR%.1f|sig%.3f|rep%d|dt%.4g]', ...
                        s.M,s.SNRdB,s.sigmaX2,iR,s.selectedDt);
                    caseSpecs{q}=s;
                end
            end
        end
    end

    taskCase=zeros(nTasks,1); taskRestart=zeros(nTasks,1); taskPriority=zeros(nTasks,1);
    checkpointPaths=cell(nTasks,1); identities=cell(nTasks,1); q=0;
    for c=1:nCases
        s=caseSpecs{c};
        for r=1:s.nStarts
            q=q+1; taskCase(q)=c; taskRestart(q)=r;
            taskPriority(q)=1e6*s.M + 1e4*(o.loghDt/max(s.selectedDt,eps)) + ...
                100*s.sigmaX2 + s.SNRdB;
            identities{q}=struct('M',s.M,'SNRdB',s.SNRdB,'SigmaX2',s.sigmaX2, ...
                'Replicate',s.iR,'Restart',r,'Seed',double(s.seed));
            if o.Checkpoint
                checkpointPaths{q}=fullfile(checkpointDir,checkpoint_filename(s,r));
            else
                checkpointPaths{q}='';
            end
        end
    end
    if q~=nTasks,error('sim_revision_AMI_vs_SNR_checkpointed_flat:TaskCount','Internal task-count mismatch.');end

    originalTaskIndex=(1:nTasks).';
    if o.LongestFirst
        [~,order]=sortrows([-taskPriority originalTaskIndex],[1 2]);
    else
        order=originalTaskIndex;
    end

    doPar=logical(o.UseParallel);
    if doPar
        try, doPar=license('test','Distrib_Computing_Toolbox')~=0; catch, doPar=false; end
    end
    pool=[];
    if doPar
        pool=gcp('nocreate');
        if isempty(pool)
            if isempty(o.ParallelWorkers),pool=parpool('Processes');
            else,pool=parpool('Processes',double(o.ParallelWorkers));end
        elseif ~isempty(o.ParallelWorkers) && pool.NumWorkers~=double(o.ParallelWorkers)
            fprintf('Existing pool has %d workers; requested %d. Reusing existing pool.\n', ...
                pool.NumWorkers,double(o.ParallelWorkers));
        end
        fprintf('Flat restart parallelism: %d process workers\n',pool.NumWorkers);
    else
        fprintf('Flat restart parallelism: disabled\n');
    end

    flatResults=cell(nTasks,1); schedulerTimer=tic;
    if doPar
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~)tick());
        parfor jj=1:nTasks
            t=order(jj); c=taskCase(t); r=taskRestart(t); spec=caseSpecs{c};
            cfg=build_case_cfg(spec,o,saProfile);
            flatResults{jj}=sa_restart_task_checkpointed(cfg,spec.baseline.x,spec.seed,r, ...
                spec.label,checkpointPaths{t},signatureText,identities{t});
            send(Q,1);
        end
    else
        for jj=1:nTasks
            t=order(jj); c=taskCase(t); r=taskRestart(t); spec=caseSpecs{c};
            cfg=build_case_cfg(spec,o,saProfile);
            flatResults{jj}=sa_restart_task_checkpointed(cfg,spec.baseline.x,spec.seed,r, ...
                spec.label,checkpointPaths{t},signatureText,identities{t});
            fprintf('  completed restart task %d/%d\n',jj,nTasks);
        end
    end
    schedulerWallSeconds=toc(schedulerTimer);

    reused=false(nTasks,1);
    for jj=1:nTasks,reused(jj)=logical(flatResults{jj}.checkpointReused);end
    nReused=sum(reused); nExecuted=nTasks-nReused;

    restartGrid=cell(nCases,1);
    for c=1:nCases,restartGrid{c}=cell(1,caseSpecs{c}.nStarts);end
    for jj=1:nTasks
        t=order(jj); c=taskCase(t); r=taskRestart(t);
        rc=restartGrid{c}; rc{r}=flatResults{jj}; restartGrid{c}=rc;
    end

    caseResults=cell(nCases,1); caseFiles=cell(nCases,1);
    for c=1:nCases
        rc=restartGrid{c}; starts=vertcat(rc{:}); spec=caseSpecs{c};
        [saOut,starts]=sa_finalize_restart_results(starts,spec.label);
        [caseResults{c},caseFiles{c}]=assemble_case(spec,starts,saOut,o,saProfile,runMeta,caseDir);
    end

    caseGrid=cell(nM,nSig,nSNR,nRep); fileGrid=cell(nM,nSig,nSNR,nRep);
    for c=1:nCases
        s=caseSpecs{c};
        caseGrid{s.iM,s.iSig,s.iSNR,s.iR}=caseResults{c};
        fileGrid{s.iM,s.iSig,s.iSNR,s.iR}=caseFiles{c};
    end
    summary=build_summary(caseGrid,baseline,MVec,sigVec,snrVec,nRep);

    ruleVersion='not-applicable';
    for k=1:numel(dtMeta)
        if isfield(dtMeta{k},'ruleVersion'),ruleVersion=dtMeta{k}.ruleVersion;break;end
    end

    bundle=struct();
    bundle.meta=runMeta; bundle.experiment='revised_AMI_vs_SNR_checkpointed_flat';
    bundle.run=struct('id',runId,'label',char(string(o.RunLabel)),'outputDirectory',outDir);
    bundle.axes=struct('M',MVec,'P_avg',Pavg,'SNR_dB',snrVec,'sigma_X_sq',sigVec,'replicate',1:nRep);
    bundle.method=struct('minGap',double(o.minGap),'pinZero',logical(o.pinZero), ...
        'fastFadingMethod',char(fastMethod),'loghDtPolicy',char(dtPolicy), ...
        'loghDtBase',double(o.loghDt),'loghDtRefined',double(o.loghDtRefined), ...
        'loghDtRuleVersion',ruleVersion,'loghDtGrid',dtGrid, ...
        'nRefinedPhysicalCases',sum(refinedMask(:)),'loghSpanSigma',double(o.loghSpanSigma), ...
        'loghYBlockSize',double(o.loghYBlockSize),'yGridMode',char(yMode),'yGridScale',double(o.YGridScale));
    bundle.sa=struct('maxIter',double(o.saMaxIter),'restartPolicy',char(restartPolicy), ...
        'startsByM',startsByM,'restartPolicyVersion',restart_policy_version(restartMeta), ...
        'nReplicates',nRep,'T0',saProfile.T0,'Tf',saProfile.Tf, ...
        'baseStd0',saProfile.baseStd0,'itersPerTemp',saProfile.itersPerTemp, ...
        'baseSeed',double(o.baseSeed),'historyEvery',double(o.historyEvery));
    bundle.parallel=struct('granularity','restart-task','useParallel',doPar, ...
        'workers',worker_count(pool,doPar),'longestFirst',logical(o.LongestFirst), ...
        'nRestartTasks',nTasks,'schedulerWallSeconds',schedulerWallSeconds,'noNestedParallelism',true);
    bundle.checkpoint=struct('enabled',logical(o.Checkpoint),'resumeRequested',resumeRequested, ...
        'manifestPath',manifestPath,'checkpointDirectory',checkpointDir, ...
        'baselineReused',baselineReused,'nReused',nReused,'nExecuted',nExecuted, ...
        'signatureText',signatureText);
    bundle.baseline=baseline; bundle.caseFiles=fileGrid; bundle.summary=summary;

    matPath=fullfile(outDir,'revisedAMIvsSNR_checkpointed_flat.mat');
    save(matPath,'bundle','-v7.3');
    if o.WriteCSV
        writetable(summary.caseTable,fullfile(outDir,'checkpointed_flat_case_summary.csv'));
    end

    fprintf('\nSaved I.11-ready bundle: %s\n',matPath);
    fprintf('Scheduler wall time: %.2f h | tasks=%d | workers=%d\n', ...
        schedulerWallSeconds/3600,nTasks,bundle.parallel.workers);
    fprintf('Checkpoint reuse: %d reused | %d executed | baseline reused=%d\n', ...
        nReused,nExecuted,baselineReused);
    fprintf('Max fast-validator gap: %.3e | selection changes: %d\n', ...
        max(summary.caseTable.MaxAbsFastValidationGap,[],'omitnan'),sum(summary.caseTable.SelectionChanged));

    function tick()
        doneCount=doneCount+1;
        fprintf('  completed restart task %d/%d\n',doneCount,nTasks);
    end
end

function signature=build_signature(MVec,snrVec,sigVec,Pavg,startsByM,o)
    signature=struct();
    signature.schedulerVersion='I11-checkpointed-flat-v1';
    signature.MVec=MVec; signature.SNRVec=snrVec; signature.TurbulenceVec=sigVec;
    signature.P_avg=Pavg; signature.minGap=double(o.minGap); signature.pinZero=logical(o.pinZero);
    signature.FastFadingMethod=char(string(o.FastFadingMethod)); signature.ghN_h=double(o.ghN_h);
    signature.LoghDtPolicy=char(string(o.LoghDtPolicy)); signature.loghDt=double(o.loghDt);
    signature.loghDtRefined=double(o.loghDtRefined); signature.loghSpanSigma=double(o.loghSpanSigma);
    signature.loghYBlockSize=double(o.loghYBlockSize); signature.YGridMode=char(string(o.YGridMode));
    signature.YGridScale=double(o.YGridScale); signature.dtRuleVersion='I10.3-v1';
    signature.runtimeTuningVersion='I10.2.2-v1';
    signature.saMaxIter=double(o.saMaxIter); signature.RestartPolicy=char(string(o.RestartPolicy));
    signature.saNStartsFixed=double(o.saNStarts); signature.startsByM=startsByM;
    signature.nReplicates=double(o.nReplicates); signature.T0=double(o.T0); signature.Tf=double(o.Tf);
    signature.BaseStd0=double(o.BaseStd0); signature.ItersPerTemp=double(o.ItersPerTemp);
    signature.baseSeed=double(o.baseSeed); signature.historyEvery=double(o.historyEvery);
end

function [dtGrid,dtMeta]=resolve_dt_grid(MVec,sigVec,snrVec,o,fastMethod)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec);
    dtGrid=repmat(double(o.loghDt),nM,nSig,nSNR); dtMeta=cell(nM,nSig,nSNR);
    for iM=1:nM
        for iSig=1:nSig
            for iSNR=1:nSNR
                if fastMethod=="logh"
                    [dtGrid(iM,iSig,iSNR),dtMeta{iM,iSig,iSNR}]=select_logh_dt( ...
                        MVec(iM),snrVec(iSNR),sigVec(iSig),'Policy',o.LoghDtPolicy, ...
                        'BaseDt',o.loghDt,'RefinedDt',o.loghDtRefined,'FixedDt',o.loghDt);
                else
                    dtMeta{iM,iSig,iSNR}=struct('policy','not-applicable', ...
                        'ruleVersion','not-applicable','selectedDt',NaN,'isRefined',false,'reason','GH evaluator');
                end
            end
        end
    end
end

function baseline=build_baselines(MVec,sigVec,snrVec,Pavg,o,xMaxBounds,dtGrid,dtMeta,saProfile,startsByM)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec); baseline=cell(nM,nSig,nSNR);
    fprintf('Validating deterministic Uniform-PAM baselines (%d physical points)...\n',nM*nSig*nSNR);
    for iM=1:nM
        x=define_constellation(MVec(iM),Pavg,true,"mean");
        for iSig=1:nSig
            for iSNR=1:nSNR
                s=struct('M',MVec(iM),'SNRdB',snrVec(iSNR),'sigmaX2',sigVec(iSig), ...
                    'selectedDt',dtGrid(iM,iSig,iSNR),'xMaxBound',xMaxBounds(iM), ...
                    'seed',1,'nStarts',startsByM(iM));
                cfg=build_case_cfg(s,o,saProfile);
                AMI_functions.assert_constellation_feasible(x,cfg,1e-10);
                t=tic; val=cfg.AMI_Validator(x(:)); valTime=toc(t); fast=cfg.AMI_Evaluator(x(:));
                m=dtMeta{iM,iSig,iSNR};
                baseline{iM,iSig,iSNR}=struct('x',x(:),'amiFast',fast,'amiValidated',val, ...
                    'fastValidationGap',fast-val,'validationRuntime',valTime, ...
                    'loghDtSelected',dtGrid(iM,iSig,iSNR),'loghDtPolicy',m.policy, ...
                    'loghDtReason',m.reason,'loghDtIsRefined',m.isRefined);
            end
        end
    end
end

function cfg=build_case_cfg(s,o,saProfile)
    cfg=build_fso_config(s.M,o.P_avg,s.SNRdB,s.sigmaX2, ...
        'FastFadingMethod',o.FastFadingMethod,'ghN_h',o.ghN_h, ...
        'loghDt',s.selectedDt,'loghSpanSigma',o.loghSpanSigma, ...
        'loghYBlockSize',o.loghYBlockSize,'YGridMode',o.YGridMode,'YGridScale',o.YGridScale, ...
        'xMaxBound',s.xMaxBound,'saMaxIter',o.saMaxIter,'saNStarts',s.nStarts, ...
        'saUseParallel',false,'minGap',o.minGap,'pinZero',o.pinZero, ...
        'seedInit',s.seed,'logEvery',0,'historyEvery',o.historyEvery);
    cfg.SA.T0=saProfile.T0; cfg.SA.Tf=saProfile.Tf; cfg.SA.baseStd0=saProfile.baseStd0;
    cfg.SA.itersPerTemp=saProfile.itersPerTemp;
    cfg.SA.nBlocks=ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate=exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);
end

function [c,path]=assemble_case(s,starts,saOut,o,saProfile,runMeta,caseDir)
    win=starts(saOut.bestStartValidated);
    k95=convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.95);
    k99=convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.99);
    vals=[starts.bestMIValidated]; totalFast=numel(starts)*double(o.saMaxIter);
    optimizerSaved=saOut; if isfield(optimizerSaved,'bestRun'),optimizerSaved=rmfield(optimizerSaved,'bestRun');end
    cfg=build_case_cfg(s,o,saProfile);
    taskSeconds=[starts.taskWallSeconds]; reused=sum([starts.checkpointReused]);

    c=struct(); c.meta=runMeta;
    c.case=struct('M',s.M,'P_avg',double(o.P_avg),'SNR_dB',s.SNRdB, ...
        'sigma_X_sq',s.sigmaX2,'replicate',s.iR,'seed',s.seed);
    c.integration=struct('loghDtSelected',s.selectedDt,'loghDtPolicy',s.dtMeta.policy, ...
        'loghDtReason',s.dtMeta.reason,'loghDtIsRefined',s.dtMeta.isRefined, ...
        'loghDtRuleVersion',s.dtMeta.ruleVersion);
    c.sa=struct('T0',cfg.SA.T0,'Tf',cfg.SA.Tf,'baseStd0',cfg.SA.baseStd0, ...
        'itersPerTemp',cfg.SA.itersPerTemp,'coolingRate',cfg.SA.coolingRate, ...
        'maxIter',cfg.SA.maxIter,'nStarts',numel(starts),'restartPolicyVersion',s.restartMeta.ruleVersion);
    c.config=fso_result_utils.config_snapshot(cfg); c.baseline=s.baseline; c.starts=starts; c.optimizer=optimizerSaved;
    c.winner=struct('startIndex',saOut.bestStartValidated,'x',saOut.bestXValidated(:), ...
        'amiFast',saOut.validatedWinnerFastMI,'amiValidated',saOut.bestMIValidated, ...
        'gainBits',saOut.bestMIValidated-s.baseline.amiValidated,'acceptanceRate',win.acceptanceRate, ...
        'convergenceIter95',k95,'convergenceIter99',k99,'history',win.history);
    c.repeatability=struct('restartValidatedMean',mean(vals,'omitnan'), ...
        'restartValidatedStd',std(vals,0,'omitnan'),'restartValidatedMin',min(vals,[],'omitnan'), ...
        'restartValidatedMax',max(vals,[],'omitnan'),'restartValidatedRange',max(vals)-min(vals));
    c.validation=struct('selectionChanged',saOut.selectionChanged,'fastWinnerIndex',saOut.bestStartFast, ...
        'validatedWinnerIndex',saOut.bestStartValidated, ...
        'validatedAdvantageOverFastWinner',saOut.validatedAdvantageOverFastWinner, ...
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap,'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    c.geometry=struct('actualMinGap',safe_min_diff(saOut.bestXValidated), ...
        'xMax',max(saOut.bestXValidated),'meanPower',mean(saOut.bestXValidated));
    c.runtime=struct('restartTaskMaxSeconds',max(taskSeconds),'restartTaskSumSeconds',sum(taskSeconds), ...
        'saSecondsSummed',saOut.totalSARuntime,'validationSeconds',saOut.validationTotalRuntime, ...
        'fastEvaluations',totalFast,'fastEvaluationsPerSecond',totalFast/max(saOut.totalSARuntime,eps), ...
        'checkpointReusedRestarts',reused);

    fname=sprintf('case_M%d_SNR%s_sig%s_rep%02d.mat',s.M,num_token(s.SNRdB,2),num_token(s.sigmaX2,4),s.iR);
    path=fullfile(caseDir,fname); save(path,'c','-v7.3');
end

function S=build_summary(caseGrid,baseline,MVec,sigVec,snrVec,nRep)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec);
    sz4=[nM nSig nSNR nRep]; sz3=[nM nSig nSNR];
    names={'amiGSValidated','amiGSFast','gainBits','convergenceIter99','restartValidatedStd', ...
        'maxAbsFastValidationGap','selectionChanged','restartTaskMaxSeconds'};
    S.raw=struct(); for k=1:numel(names),S.raw.(names{k})=nan(sz4);end
    S.xWinner=cell(sz4); S.amiPAMValidated=nan(sz3); S.loghDtSelected=nan(sz3);
    rows=cell(prod(sz4),1); q=0;
    for iM=1:nM
        for iSig=1:nSig
            for iSNR=1:nSNR
                b=baseline{iM,iSig,iSNR}; S.amiPAMValidated(iM,iSig,iSNR)=b.amiValidated;
                S.loghDtSelected(iM,iSig,iSNR)=b.loghDtSelected;
                for r=1:nRep
                    c=caseGrid{iM,iSig,iSNR,r};
                    S.raw.amiGSValidated(iM,iSig,iSNR,r)=c.winner.amiValidated;
                    S.raw.amiGSFast(iM,iSig,iSNR,r)=c.winner.amiFast;
                    S.raw.gainBits(iM,iSig,iSNR,r)=c.winner.gainBits;
                    S.raw.convergenceIter99(iM,iSig,iSNR,r)=c.winner.convergenceIter99;
                    S.raw.restartValidatedStd(iM,iSig,iSNR,r)=c.repeatability.restartValidatedStd;
                    S.raw.maxAbsFastValidationGap(iM,iSig,iSNR,r)=c.validation.maxAbsFastValidationGap;
                    S.raw.selectionChanged(iM,iSig,iSNR,r)=double(c.validation.selectionChanged);
                    S.raw.restartTaskMaxSeconds(iM,iSig,iSNR,r)=c.runtime.restartTaskMaxSeconds;
                    S.xWinner{iM,iSig,iSNR,r}=c.winner.x;
                    q=q+1;
                    rows{q}=struct('M',MVec(iM),'SigmaX2',sigVec(iSig),'SNRdB',snrVec(iSNR), ...
                        'Replicate',r,'NStarts',c.sa.nStarts,'Seed',double(c.case.seed), ...
                        'SelectedDt',c.integration.loghDtSelected,'AMIValidated',c.winner.amiValidated, ...
                        'AMIFast',c.winner.amiFast,'GainBits',c.winner.gainBits,'K99',c.winner.convergenceIter99, ...
                        'RestartValidatedStd',c.repeatability.restartValidatedStd, ...
                        'MaxAbsFastValidationGap',c.validation.maxAbsFastValidationGap, ...
                        'SelectionChanged',logical(c.validation.selectionChanged), ...
                        'RestartTaskMaxSeconds',c.runtime.restartTaskMaxSeconds, ...
                        'CheckpointReusedRestarts',c.runtime.checkpointReusedRestarts);
                end
            end
        end
    end
    rows=rows(1:q); S.caseTable=struct2table(vertcat(rows{:}));
    S.mean=struct(); S.stdAcrossReplicates=struct(); f=fieldnames(S.raw);
    for k=1:numel(f)
        S.mean.(f{k})=mean(S.raw.(f{k}),4,'omitnan');
        S.stdAcrossReplicates.(f{k})=std(S.raw.(f{k}),0,4,'omitnan');
    end
end

function name=checkpoint_filename(s,r)
    name=sprintf('restart_M%d_SNR%s_sig%s_rep%02d_r%02d.mat', ...
        s.M,num_token(s.SNRdB,2),num_token(s.sigmaX2,4),s.iR,r);
end
function v=restart_policy_version(meta)
    if isempty(meta),v='unknown';return;end
    vals=cellfun(@(m) string(m.ruleVersion),meta,'UniformOutput',false);
    vals=unique(string([vals{:}]));
    if numel(vals)==1,v=char(vals);else,v=strjoin(cellstr(vals),',');end
end
function k=convergence_iteration(history,initialMI,bestMI,fraction)
    improvement=bestMI-initialMI;
    if ~(isfinite(improvement)&&improvement>0),k=0;return;end
    target=initialMI+fraction*improvement; idx=find(history.bestMI>=target,1,'first');
    if isempty(idx),k=history.iter(end);else,k=history.iter(idx);end
end
function validate_spacing_feasibility(MVec,Pavg,minGap)
    req=(MVec-1)*minGap/2; bad=find(req>Pavg+100*eps(max(1,Pavg)),1);
    if ~isempty(bad),error('sim_revision_AMI_vs_SNR_checkpointed_flat:InfeasibleGap', ...
        'M=%d and d_min=%.6g are infeasible for P_avg=%.6g.',MVec(bad),minGap,Pavg);end
end
function n=worker_count(pool,doPar)
    if doPar && ~isempty(pool),n=pool.NumWorkers;else,n=1;end
end
function d=safe_min_diff(x)
    x=sort(real(x(:))); if numel(x)<2,d=NaN;else,d=min(diff(x));end
end
function cleanup_tmp(path)
    if exist(path,'file')==2,try,delete(path);catch,end,end
end
function token=sanitize_label(v)
    token=strtrim(char(string(v))); if isempty(token),token='run';end
    token=regexprep(token,'[^A-Za-z0-9._-]+','_');
end
function [d,id]=make_unique_dir(parent,id)
    if ~exist(parent,'dir'),mkdir(parent);end
    d=fullfile(parent,id); k=1;
    while exist(d,'dir'),k=k+1;d=fullfile(parent,sprintf('%s_%02d',id,k));end
    [ok,msg]=mkdir(d);if ~ok,error('sim_revision_AMI_vs_SNR_checkpointed_flat:mkdir','%s',msg);end
    [~,id]=fileparts(d);
end
function token=num_token(v,n)
    token=sprintf(['%0.' num2str(n) 'f'],double(v));token=strrep(token,'-','m');token=strrep(token,'.','p');
end
function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=nonnegative_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;end
function tf=positive_integer(v),tf=positive_scalar(v)&&mod(v,1)==0;end
function tf=numeric_vector(v),tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v));end
function tf=nonnegative_vector(v),tf=numeric_vector(v)&&all(v>=0);end
function tf=valid_mvec(v),tf=numeric_vector(v)&&all(v>=2)&&all(mod(v,1)==0);end
