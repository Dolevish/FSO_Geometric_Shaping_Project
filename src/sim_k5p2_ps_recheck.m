function O = sim_k5p2_ps_recheck(B,original,casePlan,taskPlan,F,prodPath,i13dPath,outDir,varargin)
%SIM_K5P2_PS_RECHECK Optional independent 0.3 PS reproduction.
% Invoked ONLY when make_k5p2_fig7('RerunPatternSearch',true).
% The original I.13D / K5 values and publication figure remain unchanged.
% Same worker, starts, constraints, validator, and 4001 evaluation budget.
% Resume is automatic for validated per-worker checkpoints and a matched
% scientific signature. Dynamic 8-process dispatch; no nested parallelism.

    p=inputParser;
    addParameter(p,'ParallelWorkers',8, ...
        @(n)isnumeric(n)&&isscalar(n)&&isfinite(n)&&n>=1&&n==floor(n));
    parse(p,varargin{:}); o=p.Results;
    assert(isequal(F.MVec,8) && isequal(F.SNRVec,20) && ...
        isequal(F.TurbulenceVec,0.3) && ...
        F.P_avg==1 && F.minGap==0.01 && F.pinZero && ...
        F.RestartsByM==6 && F.MaxIterByM==4000 && ...
        F.MaxEvalsPerStart==4001 && height(casePlan)==2 && ...
        height(taskPlan)==12,'K5p2:RecheckPolicy','Frozen policy changed.');

    subDir=fullfile(outDir,'ps_recheck');
    cpDir=fullfile(subDir,'checkpoints');
    if exist(cpDir,'dir')~=7,mkdir(cpDir);end
    P=i11p1_production_config();
    stage=B.stages.standard;
    iM=find(stage.axes.M==8,1);
    iT=find(abs(stage.axes.sigma_X_sq-.3)<1e-12,1);
    iS=find(abs(stage.axes.SNR_dB-20)<1e-12,1);
    base=stage.baseline{iM,iT,iS};
    specs=cell(2,1);
    origCases=sortrows(original.caseSummary( ...
        original.caseSummary.M==8 & original.caseSummary.SNRdB==20 & ...
        abs(original.caseSummary.SigmaR2-.3)<1e-12,:), 'Replicate');
    for c=1:2
        C=casePlan(c,:);
        if abs(base.amiValidated-C.PAMValidated)>1e-9 || ...
            C.Seed~=origCases.Seed(c) || ...
            C.MaxEvalsPerStart~=4001 || C.SelectedDt~=0.01
            error('K5p2:Baseline','Baseline, seed or evaluator discrepancy.');
        end
        cfg=build_fso_config(8,F.P_avg,20,0.3, ...
            'FastFadingMethod',P.FastFadingMethod, ...
            'loghDt',C.SelectedDt, ...
            'loghSpanSigma',P.loghSpanSigma, ...
            'loghYBlockSize',P.loghYBlockSize, ...
            'YGridMode',P.YGridMode,'YGridScale',P.YGridScale, ...
            'minGap',F.minGap,'pinZero',F.pinZero, ...
            'saMaxIter',C.MaxIterSA,'saNStarts',C.NStarts, ...
            'seedInit',C.Seed,'logEvery',0, ...
            'historyEvery',P.historyEvery);
        AMI_functions.assert_constellation_feasible(base.x(:),cfg,1e-10);
        specs{c}=struct('cfg',cfg,'xBaseline',base.x(:), ...
            'seed',uint32(C.Seed),'maxEvals',C.MaxEvalsPerStart);
    end
    signature=struct('version','K5.2-PS-recheck-v1', ...
        'production',canonical(prodPath),'originalBenchmark',canonical(i13dPath), ...
        'sciencePolicy',P,'optimizerPolicy',F.PatternSearch, ...
        'point',[8 20 0.3],'repSeeds',casePlan.Seed.', ...
        'dt',casePlan.SelectedDt.','restarts',6,'maxEvals',4001);
    sig=char(jsonencode(signature));
    manifest=fullfile(subDir,'K5p2_PS_recheck_manifest.mat');
    if exist(manifest,'file')==2
        X=load(manifest,'sig');
        if ~isfield(X,'sig')||~strcmp(X.sig,sig)
            error('K5p2:RecheckResume','Scientific signature mismatch.');
        end
    else
        save(manifest,'sig','signature');
    end
    tasks=cell(height(taskPlan),1);
    for j=1:numel(tasks)
        [c,t,id,cp]=taskargs(j,casePlan,taskPlan,cpDir);
        if exist(cp,'file')~=2,continue;end
        X=load(cp,'checkpoint');
        if ~isfield(X,'checkpoint') || ...
                ~strcmp(X.checkpoint.signatureText,sig) || ...
                ~isequaln(X.checkpoint.identity,id)
            error('K5p2:RecheckCheckpoint','Wrong checkpoint: %s',cp);
        end
        validate_result(X.checkpoint.result,specs{c}.cfg,id);
        tasks{j}=X.checkpoint.result;
    end
    pending=find(cellfun(@isempty,tasks));
    fprintf('K5.2 optional PS recheck | reused=%d/12 | pending=%d\n', ...
        numel(tasks)-numel(pending),numel(pending));
    if ~isempty(pending)
        if o.ParallelWorkers==1
            for k=1:numel(pending)
                j=pending(k);
                [c,t,id,cp]=taskargs(j,casePlan,taskPlan,cpDir);
                fprintf('K5.2 PS START %d/%d | %s\n',k,numel(pending),id.IdentityKey);
                s=specs{c};
                r=i13d_patternsearch_restart_task(s.cfg,s.xBaseline, ...
                    s.seed,t.Restart,s.maxEvals,F.PatternSearch,id,cp,sig);
                validate_result(r,s.cfg,id);
                tasks{j}=r;
                fprintf('K5.2 PS DONE %d/%d | %s | %.2fs | val=%.8f\n', ...
                    k,numel(pending),id.IdentityKey,r.taskWallSeconds,r.bestMIValidated);
            end
        else
            pp=gcp('nocreate');
            if ~isempty(pp),delete(pp);end
            pp=parpool('Processes',double(o.ParallelWorkers));
            guard=onCleanup(@() delete_pool(pp)); %#ok<NASGU>
            futures=parallel.FevalFuture.empty(0,1);
            indices=zeros(0,1);
            next=1;done=0;
            try
                while numel(futures)<pp.NumWorkers && next<=numel(pending)
                    j=pending(next);
                    futures(end+1,1)=submit(pp,j,casePlan,taskPlan,cpDir,specs,F,sig); %#ok<AGROW>
                    indices(end+1)=j; %#ok<AGROW>
                    fprintf('K5.2 PS START %d/%d | %s | active=%d/%d\n', ...
                        next,numel(pending),taskPlan.IdentityKey(j),numel(futures),pp.NumWorkers);
                    next=next+1;
                end
                while ~isempty(futures)
                    [pos,r]=fetchNext(futures);
                    j=indices(pos);
                    futures(pos)=[];indices(pos)=[];
                    [c,~,id,~]=taskargs(j,casePlan,taskPlan,cpDir);
                    validate_result(r,specs{c}.cfg,id);
                    tasks{j}=r;done=done+1;
                    fprintf('K5.2 PS DONE %d/%d | %s | %.2fs | val=%.8f\n', ...
                        done,numel(pending),id.IdentityKey,r.taskWallSeconds,r.bestMIValidated);
                    if next<=numel(pending)
                        j=pending(next);
                        futures(end+1,1)=submit(pp,j,casePlan,taskPlan,cpDir,specs,F,sig); %#ok<AGROW>
                        indices(end+1)=j; %#ok<AGROW>
                        fprintf('K5.2 PS START %d/%d | %s | active=%d/%d\n', ...
                            next,numel(pending),taskPlan.IdentityKey(j),numel(futures),pp.NumWorkers);
                        next=next+1;
                    end
                end
            catch ME
                if ~isempty(futures),try,cancel(futures);catch,end,end
                rethrow(ME);
            end
        end
    end

    Replicate=zeros(12,1);Restart=zeros(12,1);
    PSFast=zeros(12,1);PSValidated=zeros(12,1);
    FunctionEvaluations=zeros(12,1);TaskWallSeconds=zeros(12,1);
    for j=1:12
        [c,~,id,~]=taskargs(j,casePlan,taskPlan,cpDir);
        r=tasks{j};validate_result(r,specs{c}.cfg,id);
        Replicate(j)=r.Replicate;Restart(j)=r.Restart;
        PSFast(j)=r.bestMIFast;PSValidated(j)=r.bestMIValidated;
        FunctionEvaluations(j)=r.functionEvaluations;
        TaskWallSeconds(j)=r.taskWallSeconds;
    end
    T=sortrows(table(Replicate,Restart,PSFast,PSValidated, ...
        FunctionEvaluations,TaskWallSeconds),{'Replicate','Restart'});
    winner=zeros(2,1);
    for rep=1:2
        x=T.PSValidated(T.Replicate==rep);
        if numel(x)~=6,error('K5p2:RecheckCount','Missing replicate restart.');end
        winner(rep)=max(x);
    end
    old=origCases.PSValidatedCase;
    perRep=table([1;2],old,winner,winner-old, ...
        'VariableNames',{'Replicate','FrozenI13D','IndependentRecheck','Delta'});
    writetable(T,fullfile(subDir,'K5p2_rerun_patternsearch_restarts.csv'));
    writetable(perRep,fullfile(subDir,'K5p2_rerun_vs_I13D_replicates.csv'));
    summary=table(mean(old),mean(winner),mean(winner)-mean(old), ...
        std(old,0),std(winner,0), ...
        'VariableNames',{'FrozenPSMean','RecheckPSMean','DeltaMean', ...
        'FrozenPSStd','RecheckPSStd'});
    writetable(summary,fullfile(subDir,'K5p2_rerun_vs_I13D_summary.csv'));
    O=struct('version','K5.2-PS-recheck-v1','sourceOriginalI13D',i13dPath, ...
        'sourceProduction',prodPath,'outputDirectory',subDir, ...
        'restartResults',T,'comparison',perRep,'summary',summary, ...
        'noOverwriteOfPublicationFigure',true);
    save(fullfile(subDir,'K5p2_rerun_results.mat'),'O','-v7.3');
    fprintf('K5.2 PS recheck: frozen=%.8f, new=%.8f, delta=%+.3e\n', ...
        mean(old),mean(winner),mean(winner)-mean(old));
end

function f=submit(pp,j,casePlan,taskPlan,cpDir,specs,F,sig)
    [c,t,id,cp]=taskargs(j,casePlan,taskPlan,cpDir);
    s=specs{c};
    f=parfeval(pp,@i13d_patternsearch_restart_task,1,s.cfg,s.xBaseline, ...
        s.seed,t.Restart,s.maxEvals,F.PatternSearch,id,cp,sig);
end
function [c,t,id,cp]=taskargs(j,casePlan,taskPlan,cpDir)
    t=taskPlan(j,:);c=t.TaskCaseIndex;
    x=casePlan(c,:);
    id=struct('IdentityKey',char(t.IdentityKey),'M',x.M, ...
        'SNRdB',x.SNRdB,'SigmaR2',x.SigmaR2, ...
        'Replicate',x.Replicate,'Restart',t.Restart,'Seed',double(x.Seed));
    cp=fullfile(cpDir,[id.IdentityKey '.mat']);
end
function validate_result(r,cfg,id)
    if isempty(r) || ~isstruct(r) || ~isfield(r,'success') || ~r.success
        if isstruct(r)&&isfield(r,'errorReport')
            error('K5p2:RecheckWorker','%s',char(r.errorReport));
        end
        error('K5p2:RecheckWorker','Restart failed.');
    end
    if string(r.IdentityKey)~=string(id.IdentityKey) || ...
            r.M~=8 || r.SNRdB~=20 || abs(r.SigmaR2-.3)>1e-12 || ...
            r.Replicate~=id.Replicate || r.Restart~=id.Restart || ...
            r.masterSeed~=id.Seed || ...
            r.functionEvaluations>4001 || r.functionEvaluations<1 || ...
            ~isfinite(r.bestMIValidated) || ...
            ~isfinite(r.bestMIFast) || ...
            abs(r.bestMIValidated-r.bestMIFast)>1e-3 || ...
            r.postHocFastAuditEvaluations~=1
        error('K5p2:RecheckInvalid','Invalid worker identity, budget or validation.');
    end
    AMI_functions.assert_constellation_feasible(r.xBest(:),cfg,1e-10);
end
function delete_pool(pp)
    if ~isempty(pp),try,delete(pp);catch,end,end
end
function s=canonical(path)
    s=char(java.io.File(char(string(path))).getCanonicalPath());
end
