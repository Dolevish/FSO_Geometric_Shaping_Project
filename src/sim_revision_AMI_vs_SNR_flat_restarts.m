function bundle = sim_revision_AMI_vs_SNR_flat_restarts(varargin)
%SIM_REVISION_AMI_VS_SNR_FLAT_RESTARTS  Flat restart-task production scheduler.
%
% Commit I.10.5A removes the long-tail weakness of case-level parallelism.
% The atomic parallel job is one SA restart (including independent
% validation), identified by
%
%   (M, SNR_dB, sigma_X^2, replicate, restart).
%
% All restart jobs from all optimization cases share ONE process-pool queue.
% There is no nested parallelism. When a short case finishes, its worker can
% immediately execute a remaining restart from a heavier case such as M=32.
%
% Reproducibility:
%   - the physical-case seed is identical to sim_revision_AMI_vs_SNR;
%   - restart initializations and SA streams are reconstructed by
%     sa_restart_task using the same Threefry masterSeed/substream mapping as
%     sa_multistart;
%   - final selection remains the maximum independently VALIDATED AMI across
%     all restarts of a case.
%
% The scheduler defaults to longest-first ordering by M/refined-dt/turbulence
% to reduce the end-of-run heavy-task tail. Ordering changes only wall-clock
% scheduling; each restart's random stream is index-addressed and invariant
% to execution order.
%
% This driver is intentionally separate from the historical production
% driver until the flat scheduler passes equivalence and scaling diagnostics.

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
    addParameter(p,'WriteCSV',true,@logical_scalar);
    addParameter(p,'RunLabel','',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    parse(p,varargin{:}); o=p.Results;

    if ~(o.T0>o.Tf)
        error('sim_revision_AMI_vs_SNR_flat_restarts:TemperatureOrder','T0 must exceed Tf.');
    end
    if o.loghDtRefined>o.loghDt
        error('sim_revision_AMI_vs_SNR_flat_restarts:DtOrder','Refined dt must not exceed base dt.');
    end

    fastMethod=lower(string(o.FastFadingMethod));
    yMode=lower(string(o.YGridMode));
    dtPolicy=lower(string(o.LoghDtPolicy));
    if ~ismember(fastMethod,["logh","gh"])
        error('sim_revision_AMI_vs_SNR_flat_restarts:FastMethod','FastFadingMethod must be logh or gh.');
    end
    if ~ismember(yMode,["fixed-bound","candidate-adaptive"])
        error('sim_revision_AMI_vs_SNR_flat_restarts:YGridMode','Unsupported YGridMode.');
    end
    if ~ismember(dtPolicy,["case-adaptive","fixed"])
        error('sim_revision_AMI_vs_SNR_flat_restarts:DtPolicy','Unsupported LoghDtPolicy.');
    end
    if fastMethod=="gh" && yMode=="candidate-adaptive"
        error('sim_revision_AMI_vs_SNR_flat_restarts:AdaptiveRequiresLogh', ...
            'candidate-adaptive y-grid is supported only by the log-h evaluator.');
    end
    o.FastFadingMethod=char(fastMethod); o.YGridMode=char(yMode); o.LoghDtPolicy=char(dtPolicy);

    MVec=unique(double(o.MVec(:).'),'stable');
    snrVec=unique(double(o.SNRVec(:).'),'stable');
    sigVec=unique(double(o.TurbulenceVec(:).'),'stable');
    Pavg=double(o.P_avg); nRep=double(o.nReplicates); nStarts=double(o.saNStarts);
    validate_spacing_feasibility(MVec,Pavg,double(o.minGap));

    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec);
    nPhysical=nM*nSig*nSNR; nCases=nPhysical*nRep; nTasks=nCases*nStarts;

    xMaxBounds=nan(1,nM);
    for iM=1:nM
        xMaxBounds(iM)=fso_result_utils.feasible_xmax_bound(MVec(iM),Pavg,double(o.minGap),logical(o.pinZero));
    end

    [dtGrid,dtMeta]=resolve_dt_grid(MVec,sigVec,snrVec,o,fastMethod);
    refinedMask=cellfun(@(m) isfield(m,'isRefined')&&m.isRefined,dtMeta);

    parentDir=fullfile(char(o.ResultsRoot),'revised_ami_vs_snr_flat_restarts');
    label=sanitize_label(o.RunLabel);
    token=sprintf('I%d_S%d_R%d',o.saMaxIter,o.saNStarts,o.nReplicates);
    utc=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId]=make_unique_dir(parentDir,sprintf('%s_%s_%s',label,token,utc));
    caseDir=fullfile(outDir,'cases'); mkdir(caseDir);

    runMeta=fso_result_utils.run_metadata(mfilename);
    runMeta.runId=runId; runMeta.runLabel=char(string(o.RunLabel)); runMeta.outputDirectory=outDir;

    saProfile=struct('T0',double(o.T0),'Tf',double(o.Tf), ...
        'baseStd0',double(o.BaseStd0),'itersPerTemp',double(o.ItersPerTemp));

    fprintf('\n============================================================\n');
    fprintf('I.10.5A flat restart-task AMI-vs-SNR scheduler\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | replicates=%d\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),nRep);
    fprintf('SA=%d restarts x %d iterations | flattened tasks=%d\n',nStarts,o.saMaxIter,nTasks);
    fprintf('Scheduling: single-level restart tasks | longest-first=%d\n',logical(o.LongestFirst));
    fprintf('Refined physical points: %d/%d\n',sum(refinedMask(:)),nPhysical);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    baseline=build_baselines(MVec,sigVec,snrVec,Pavg,o,xMaxBounds,dtGrid,dtMeta,saProfile);

    caseSpecs=cell(nCases,1);
    q=0;
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
                    s.label=sprintf('[FlatRev|M%d|SNR%.1f|sig%.3f|rep%d|dt%.4g]', ...
                        s.M,s.SNRdB,s.sigmaX2,iR,s.selectedDt);
                    caseSpecs{q}=s;
                end
            end
        end
    end

    taskCase=zeros(nTasks,1); taskRestart=zeros(nTasks,1); taskPriority=zeros(nTasks,1);
    q=0;
    for c=1:nCases
        s=caseSpecs{c};
        for r=1:nStarts
            q=q+1; taskCase(q)=c; taskRestart(q)=r;
            taskPriority(q)=1e6*s.M + 1e4*(o.loghDt/max(s.selectedDt,eps)) + ...
                100*s.sigmaX2 + s.SNRdB;
        end
    end
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
            if isempty(o.ParallelWorkers), pool=parpool('Processes');
            else, pool=parpool('Processes',double(o.ParallelWorkers)); end
        end
        fprintf('Flat restart parallelism: %d process workers\n',pool.NumWorkers);
        if ~isempty(o.ParallelWorkers) && pool.NumWorkers~=double(o.ParallelWorkers)
            fprintf('  requested %d workers; existing pool has %d and is reused\n', ...
                o.ParallelWorkers,pool.NumWorkers);
        end
    else
        fprintf('Flat restart parallelism: disabled\n');
    end

    flatResults=cell(nTasks,1);
    schedulerTimer=tic;
    if doPar
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~)tick());
        parfor jj=1:nTasks
            t=order(jj); c=taskCase(t); r=taskRestart(t); spec=caseSpecs{c};
            cfg=build_case_cfg(spec,o,saProfile);
            flatResults{jj}=sa_restart_task(cfg,spec.baseline.x,spec.seed,r,spec.label);
            send(Q,1);
        end
    else
        for jj=1:nTasks
            t=order(jj); c=taskCase(t); r=taskRestart(t); spec=caseSpecs{c};
            cfg=build_case_cfg(spec,o,saProfile);
            flatResults{jj}=sa_restart_task(cfg,spec.baseline.x,spec.seed,r,spec.label);
            fprintf('  completed restart task %d/%d\n',jj,nTasks);
        end
    end
    schedulerWallSeconds=toc(schedulerTimer);

    restartGrid=cell(nCases,nStarts);
    for jj=1:nTasks
        t=order(jj);
        restartGrid{taskCase(t),taskRestart(t)}=flatResults{jj};
    end

    caseResults=cell(nCases,1); caseFiles=cell(nCases,1);
    for c=1:nCases
        starts=vertcat(restartGrid{c,:});
        spec=caseSpecs{c};
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

    bundle=struct();
    bundle.meta=runMeta; bundle.experiment='revised_AMI_vs_SNR_flat_restarts';
    bundle.run=struct('id',runId,'label',char(string(o.RunLabel)),'outputDirectory',outDir);
    bundle.axes=struct('M',MVec,'P_avg',Pavg,'SNR_dB',snrVec,'sigma_X_sq',sigVec,'replicate',1:nRep);
    ruleVersion='not-applicable';
    for k=1:numel(dtMeta)
        if isfield(dtMeta{k},'ruleVersion'), ruleVersion=dtMeta{k}.ruleVersion; break; end
    end
    bundle.method=struct('minGap',double(o.minGap),'pinZero',logical(o.pinZero), ...
        'fastFadingMethod',char(fastMethod),'loghDtPolicy',char(dtPolicy), ...
        'loghDtBase',double(o.loghDt),'loghDtRefined',double(o.loghDtRefined), ...
        'loghDtRuleVersion',ruleVersion,'loghDtGrid',dtGrid, ...
        'nRefinedPhysicalCases',sum(refinedMask(:)),'loghSpanSigma',double(o.loghSpanSigma), ...
        'yGridMode',char(yMode),'yGridScale',double(o.YGridScale));
    bundle.sa=struct('maxIter',double(o.saMaxIter),'nStarts',nStarts,'nReplicates',nRep, ...
        'T0',saProfile.T0,'Tf',saProfile.Tf,'baseStd0',saProfile.baseStd0, ...
        'itersPerTemp',saProfile.itersPerTemp,'baseSeed',double(o.baseSeed), ...
        'historyEvery',double(o.historyEvery));
    bundle.parallel=struct('granularity','restart-task','useParallel',doPar, ...
        'workers',worker_count(pool,doPar),'longestFirst',logical(o.LongestFirst), ...
        'nRestartTasks',nTasks,'schedulerWallSeconds',schedulerWallSeconds, ...
        'noNestedParallelism',true);
    bundle.baseline=baseline; bundle.caseFiles=fileGrid; bundle.summary=summary;

    matPath=fullfile(outDir,'revisedAMIvsSNR_flat_restarts.mat');
    save(matPath,'bundle','-v7.3');
    if o.WriteCSV
        writetable(summary.caseTable,fullfile(outDir,'flat_restart_case_summary.csv'));
    end

    fprintf('\nSaved flat-restart bundle: %s\n',matPath);
    fprintf('Scheduler wall time: %.2f h | tasks=%d | workers=%d\n', ...
        schedulerWallSeconds/3600,nTasks,bundle.parallel.workers);
    fprintf('Max fast-validator gap: %.3e | selection changes: %d\n', ...
        max(summary.caseTable.MaxAbsFastValidationGap,[],'omitnan'), ...
        sum(summary.caseTable.SelectionChanged));

    function tick()
        doneCount=doneCount+1;
        fprintf('  completed restart task %d/%d\n',doneCount,nTasks);
    end
end


function [dtGrid,dtMeta]=resolve_dt_grid(MVec,sigVec,snrVec,o,fastMethod)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec);
    dtGrid=repmat(double(o.loghDt),nM,nSig,nSNR); dtMeta=cell(nM,nSig,nSNR);
    for iM=1:nM
        for iSig=1:nSig
            for iSNR=1:nSNR
                if fastMethod=="logh"
                    [dtGrid(iM,iSig,iSNR),dtMeta{iM,iSig,iSNR}]=select_logh_dt( ...
                        MVec(iM),snrVec(iSNR),sigVec(iSig), ...
                        'Policy',o.LoghDtPolicy,'BaseDt',o.loghDt, ...
                        'RefinedDt',o.loghDtRefined,'FixedDt',o.loghDt);
                else
                    dtMeta{iM,iSig,iSNR}=struct('policy','not-applicable','ruleVersion','not-applicable', ...
                        'selectedDt',NaN,'isRefined',false,'reason','GH evaluator');
                end
            end
        end
    end
end


function baseline=build_baselines(MVec,sigVec,snrVec,Pavg,o,xMaxBounds,dtGrid,dtMeta,saProfile)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec); baseline=cell(nM,nSig,nSNR);
    fprintf('Validating deterministic Uniform-PAM baselines (%d physical points)...\n',nM*nSig*nSNR);
    for iM=1:nM
        x=define_constellation(MVec(iM),Pavg,true,"mean");
        for iSig=1:nSig
            for iSNR=1:nSNR
                s=struct('M',MVec(iM),'SNRdB',snrVec(iSNR),'sigmaX2',sigVec(iSig), ...
                    'selectedDt',dtGrid(iM,iSig,iSNR),'xMaxBound',xMaxBounds(iM),'seed',1);
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
        'xMaxBound',s.xMaxBound,'saMaxIter',o.saMaxIter,'saNStarts',o.saNStarts, ...
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
    vals=[starts.bestMIValidated]; totalFast=double(o.saNStarts)*double(o.saMaxIter);
    optimizerSaved=saOut; if isfield(optimizerSaved,'bestRun'), optimizerSaved=rmfield(optimizerSaved,'bestRun'); end
    cfg=build_case_cfg(s,o,saProfile);

    actualCaseSpan=max([starts.taskEndPosix])-min([starts.taskStartPosix]);
    criticalPath=max([starts.taskWallSeconds]);

    c=struct(); c.meta=runMeta;
    c.case=struct('M',s.M,'P_avg',double(o.P_avg),'SNR_dB',s.SNRdB, ...
        'sigma_X_sq',s.sigmaX2,'replicate',s.iR,'seed',s.seed);
    c.integration=struct('loghDtSelected',s.selectedDt,'loghDtPolicy',s.dtMeta.policy, ...
        'loghDtReason',s.dtMeta.reason,'loghDtIsRefined',s.dtMeta.isRefined, ...
        'loghDtRuleVersion',s.dtMeta.ruleVersion);
    c.sa=struct('T0',cfg.SA.T0,'Tf',cfg.SA.Tf,'baseStd0',cfg.SA.baseStd0, ...
        'itersPerTemp',cfg.SA.itersPerTemp,'coolingRate',cfg.SA.coolingRate, ...
        'maxIter',cfg.SA.maxIter,'nStarts',cfg.SA.nStarts);
    c.config=fso_result_utils.config_snapshot(cfg); c.baseline=s.baseline;
    c.starts=starts; c.optimizer=optimizerSaved;
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
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap, ...
        'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    c.geometry=struct('actualMinGap',safe_min_diff(saOut.bestXValidated), ...
        'xMax',max(saOut.bestXValidated),'meanPower',mean(saOut.bestXValidated));
    c.runtime=struct('wallSeconds',actualCaseSpan,'parallelCriticalPathSeconds',criticalPath, ...
        'saSecondsSummed',saOut.totalSARuntime,'validationSeconds',saOut.validationTotalRuntime, ...
        'fastEvaluations',totalFast,'fastEvaluationsPerSecond',totalFast/max(saOut.totalSARuntime,eps), ...
        'wallSecondsMeaning','span from earliest restart-task start to latest restart-task finish');

    fname=sprintf('case_M%d_SNR%s_sig%s_rep%02d.mat',s.M,num_token(s.SNRdB,2),num_token(s.sigmaX2,4),s.iR);
    path=fullfile(caseDir,fname); save(path,'c','-v7.3');
end


function S=build_summary(caseGrid,baseline,MVec,sigVec,snrVec,nRep)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec);
    sz4=[nM nSig nSNR nRep]; sz3=[nM nSig nSNR];
    names={'amiGSValidated','amiGSFast','gainBits','convergenceIter99','restartValidatedStd', ...
        'maxAbsFastValidationGap','selectionChanged','wallSeconds','parallelCriticalPathSeconds'};
    S.raw=struct(); for k=1:numel(names), S.raw.(names{k})=nan(sz4); end
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
                    S.raw.wallSeconds(iM,iSig,iSNR,r)=c.runtime.wallSeconds;
                    S.raw.parallelCriticalPathSeconds(iM,iSig,iSNR,r)=c.runtime.parallelCriticalPathSeconds;
                    S.xWinner{iM,iSig,iSNR,r}=c.winner.x;
                    q=q+1;
                    rows{q}=struct('M',MVec(iM),'SigmaX2',sigVec(iSig),'SNRdB',snrVec(iSNR), ...
                        'Replicate',r,'Seed',double(c.case.seed),'SelectedDt',c.integration.loghDtSelected, ...
                        'AMIValidated',c.winner.amiValidated,'AMIFast',c.winner.amiFast, ...
                        'GainBits',c.winner.gainBits,'K99',c.winner.convergenceIter99, ...
                        'RestartValidatedStd',c.repeatability.restartValidatedStd, ...
                        'MaxAbsFastValidationGap',c.validation.maxAbsFastValidationGap, ...
                        'SelectionChanged',logical(c.validation.selectionChanged), ...
                        'CaseWallSeconds',c.runtime.wallSeconds, ...
                        'ParallelCriticalPathSeconds',c.runtime.parallelCriticalPathSeconds);
                end
            end
        end
    end
    rows=rows(1:q); S.caseTable=struct2table(vertcat(rows{:}));
    S.mean=struct(); S.stdAcrossReplicates=struct();
    f=fieldnames(S.raw);
    for k=1:numel(f)
        S.mean.(f{k})=mean(S.raw.(f{k}),4,'omitnan');
        S.stdAcrossReplicates.(f{k})=std(S.raw.(f{k}),0,4,'omitnan');
    end
end


function k=convergence_iteration(history,initialMI,bestMI,fraction)
    improvement=bestMI-initialMI;
    if ~(isfinite(improvement)&&improvement>0), k=0; return; end
    target=initialMI+fraction*improvement;
    idx=find(history.bestMI>=target,1,'first');
    if isempty(idx), k=history.iter(end); else, k=history.iter(idx); end
end

function validate_spacing_feasibility(MVec,Pavg,minGap)
    req=(MVec-1)*minGap/2; bad=find(req>Pavg+100*eps(max(1,Pavg)),1);
    if ~isempty(bad)
        error('sim_revision_AMI_vs_SNR_flat_restarts:InfeasibleGap', ...
            'M=%d and d_min=%.6g are infeasible for P_avg=%.6g.',MVec(bad),minGap,Pavg);
    end
end

function n=worker_count(pool,doPar)
    if doPar && ~isempty(pool), n=pool.NumWorkers; else, n=1; end
end
function d=safe_min_diff(x)
    x=sort(real(x(:)));
    if numel(x)<2, d=NaN; else, d=min(diff(x)); end
end
function token=sanitize_label(v)
    token=strtrim(char(string(v)));
    if isempty(token),token='run';end
    token=regexprep(token,'[^A-Za-z0-9._-]+','_');
end
function [d,id]=make_unique_dir(parent,id)
    if ~exist(parent,'dir'),mkdir(parent);end
    d=fullfile(parent,id); k=1;
    while exist(d,'dir'),k=k+1; d=fullfile(parent,sprintf('%s_%02d',id,k));end
    [ok,msg]=mkdir(d); if ~ok,error('sim_revision_AMI_vs_SNR_flat_restarts:mkdir','%s',msg);end
    [~,id]=fileparts(d);
end
function token=num_token(v,n)
    token=sprintf(['%0.' num2str(n) 'f'],double(v)); token=strrep(token,'-','m'); token=strrep(token,'.','p');
end
function tf=text_scalar(v), tf=ischar(v)||(isstring(v)&&isscalar(v)); end
function tf=logical_scalar(v), tf=islogical(v)&&isscalar(v); end
function tf=positive_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0; end
function tf=nonnegative_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0; end
function tf=positive_integer(v), tf=positive_scalar(v)&&mod(v,1)==0; end
function tf=numeric_vector(v), tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v)); end
function tf=nonnegative_vector(v), tf=numeric_vector(v)&&all(v>=0); end
function tf=valid_mvec(v), tf=numeric_vector(v)&&all(v>=2)&&all(mod(v,1)==0); end
