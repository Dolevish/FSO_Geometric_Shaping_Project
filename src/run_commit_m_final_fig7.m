function M = run_commit_m_final_fig7(varargin)
%RUN_COMMIT_M_FINAL_FIG7  Final paired SA-vs-GPS Figure 7 experiment.
% COMMIT M: 12 distinct paired starts; classic GPSPositiveBasis2N and
% UseCompletePoll=true; all other M=8/SNR=20/sigma_R^2=0.3 settings frozen.
% Figure reports MAXIMUM independently validated AMI of the 12 runs.
% Mean/variance/std of ALL 12 outcomes are saved, not shown on the figure.
% Frozen Uniform-PAM and Zhou-GS baselines are READ, never re-optimized.
%
% Run from repository root:
%   addpath(genpath(fullfile(pwd,'src')));
%   test_commit_m_preflight();
%   plan = run_commit_m_final_fig7('DryRun',true);
%   M = run_commit_m_final_fig7();
% If frozen K5.2 CSV is not found automatically:
%   M = run_commit_m_final_fig7('ReferenceCSVPath','/exact/path/Fig7_K5p2_four_methods.csv');
% After interruption:
%   M = run_commit_m_final_fig7('ResumeDirectory','/exact/results/ieee_revision2/commit_m_final_fig7/M_...');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ReferenceCSVPath','',@text_arg);
    addParameter(p,'ResumeDirectory','',@text_arg);
    addParameter(p,'ParallelWorkers',8,@positive_int);
    addParameter(p,'UseParallel',true,@scalar_logical);
    addParameter(p,'DryRun',false,@scalar_logical);
    parse(p,varargin{:}); o=p.Results;

    P=i11p1_production_config();
    F=i13d_benchmark_config();
    F.PatternSearch.UseCompletePoll=true;
    Msize=8; snr=20; sig=0.3; nStarts=12;
    maxIter=4000; maxEvals=maxIter+1;
    assert(P.P_avg==1 && P.minGap==0.01 && P.pinZero && ...
        strcmp(P.FastFadingMethod,'logh') && strcmp(P.YGridMode,'candidate-adaptive') && ...
        P.loghDt==0.01 && P.YGridScale==1.1 && P.loghSpanSigma==8 && ...
        P.loghYBlockSize==512 && P.BaseStd0==0.15 && ...
        P.T0==0.4 && P.Tf==1e-3 && P.ItersPerTemp==50 && ...
        F.PatternSearch.UseCompletePoll && ~F.PatternSearch.UseParallel && ...
        F.MaxEvalsPerStart(F.MVec==8)==maxEvals, ...
        'CommitM:FrozenConfig','Source settings no longer match the original Figure 7.');

    cfg=build_fso_config(Msize,P.P_avg,snr,sig, ...
        'FastFadingMethod',P.FastFadingMethod, ...
        'loghDt',P.loghDt,'loghSpanSigma',P.loghSpanSigma, ...
        'loghYBlockSize',P.loghYBlockSize,'YGridMode',P.YGridMode, ...
        'YGridScale',P.YGridScale,'minGap',P.minGap,'pinZero',P.pinZero, ...
        'saMaxIter',maxIter,'saNStarts',nStarts, ...
        'seedInit',0,'logEvery',0,'historyEvery',P.historyEvery);
    cfg.SA.T0=P.T0; cfg.SA.Tf=P.Tf;
    cfg.SA.baseStd0=P.BaseStd0; cfg.SA.itersPerTemp=P.ItersPerTemp;
    cfg.SA.nBlocks=ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate=exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);
    xBaseline=define_constellation(Msize,P.P_avg,true,"mean");
    masterSeed=fso_result_utils.case_seed(P.baseSeed,Msize,snr,sig);

    Xraw=i13_shared_start_points(cfg,xBaseline,masterSeed,nStarts);
    Xfeasible=zeros(Msize,nStarts);
    for k=1:nStarts
        Xfeasible(:,k)=AMI_functions.project_constellation_1D(Xraw(:,k),cfg);
    end
    for i=1:nStarts
        for j=i+1:nStarts
            assert(norm(Xfeasible(:,i)-Xfeasible(:,j),inf)>1e-10, ...
                'CommitM:DuplicateStart','Restarts %d and %d share a projected start.',i,j);
        end
    end
    assert(norm(Xfeasible(:,1)-xBaseline(:),inf)<1e-12, ...
        'CommitM:BaselineStart','Restart 1 must be Uniform PAM.');

    signature=struct('version','CommitM-v1','M',Msize,'SNRdB',snr, ...
        'SigmaR2',sig,'Pavg',P.P_avg,'minGap',P.minGap, ...
        'pinZero',P.pinZero,'nStarts',nStarts,'masterSeed',masterSeed, ...
        'maxIterSA',maxIter,'maxEvalsGPS',maxEvals,'saT0',P.T0, ...
        'saTf',P.Tf,'saBaseStd0',P.BaseStd0, ...
        'saItersPerTemp',P.ItersPerTemp, ...
        'loghDt',P.loghDt,'loghSpanSigma',P.loghSpanSigma, ...
        'YGridMode',P.YGridMode,'YGridScale',P.YGridScale, ...
        'loghYBlockSize',P.loghYBlockSize, ...
        'PatternSearch',F.PatternSearch, ...
        'startPolicy','first=Uniform PAM; 2:12=Threefry Gaussian 0.2; paired SA/GPS');
    signatureText=char(jsonencode(signature));

    taskMethod=[repmat("SA",nStarts,1);repmat("GPS",nStarts,1)];
    taskRestart=[(1:nStarts)';(1:nStarts)'];
    taskPlan=table(taskMethod,taskRestart,'VariableNames',{'Method','Restart'});

    if o.DryRun
        M=struct('version','CommitM-v1','dryRun',true,'signature',signature, ...
            'taskPlan',taskPlan,'rawStarts',Xraw,'feasibleStarts',Xfeasible, ...
            'referenceCSV',string(o.ReferenceCSVPath), ...
            'resultsRoot',fso_result_utils.revision2_results_root());
        fprintf('Commit M dry run: 12 distinct paired starts, 24 tasks, complete GPS polling.\n');
        return;
    end

    refPath=resolve_reference_csv(o.ReferenceCSVPath);
    frozen=readtable(refPath,'TextType','string');
    names=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    assert(all(ismember({'Method','ValidatedAMI'},frozen.Properties.VariableNames)) && ...
        height(frozen)==4 && isequal(string(frozen.Method(:)),names), ...
        'CommitM:FrozenCSV','Expected the four method rows from K5.2.');
    oldValues=double(frozen.ValidatedAMI(:));
    assert(all(isfinite(oldValues)) && ...
        abs(oldValues(1)-0.9323308194)<2e-6 && ...
        abs(oldValues(2)-1.1878905650)<2e-6, ...
        'CommitM:FrozenBaselines','Frozen CSV does not match the approved M8/20/0.3 baselines.');
    baselineValidated=cfg.AMI_Validator(xBaseline(:));
    assert(abs(baselineValidated-oldValues(1))<1e-6, ...
        'CommitM:BaselineMismatch','Current validator disagrees with frozen Uniform PAM.');

    root=fullfile(fso_result_utils.revision2_results_root(),'commit_m_final_fig7');
    if strlength(string(o.ResumeDirectory))>0
        runDir=char(o.ResumeDirectory);
        assert(isfolder(runDir),'CommitM:MissingResume','ResumeDirectory is missing.');
    else
        stamp=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
        runDir=fullfile(root,['M_' stamp]);
        if ~isfolder(runDir),mkdir(runDir);end
    end
    cpDir=fullfile(runDir,'checkpoints');
    nativeDir=fullfile(runDir,'native_gps');
    if ~isfolder(cpDir),mkdir(cpDir);end
    if ~isfolder(nativeDir),mkdir(nativeDir);end

    resultFile=fullfile(runDir,'CommitM_results.mat');
    results=cell(height(taskPlan),1);
    checkpointPaths=cell(height(taskPlan),1);
    for k=1:height(taskPlan)
        key=sprintf('%s_restart_%02d',char(taskPlan.Method(k)),taskPlan.Restart(k));
        checkpointPaths{k}=fullfile(cpDir,[key '.mat']);
        if exist(checkpointPaths{k},'file')==2
            C=load(checkpointPaths{k},'checkpoint');
            assert(isfield(C,'checkpoint') && isfield(C.checkpoint,'signatureText') && ...
                strcmp(C.checkpoint.signatureText,signatureText) && ...
                strcmp(C.checkpoint.key,key) && ...
                isfield(C.checkpoint,'result'), ...
                'CommitM:CheckpointMismatch','Checkpoint identity/science mismatch: %s',key);
            results{k}=C.checkpoint.result;
            verify_run_result(results{k},taskMethod(k),taskRestart(k),Xfeasible(:,taskRestart(k)),cfg);
        end
    end
    pending=find(cellfun(@isempty,results));
    fprintf('Commit M: %d completed, %d pending | results: %s\n',numel(results)-numel(pending),numel(pending),runDir);

    if ~isempty(pending)
        if o.UseParallel
            pool=gcp('nocreate');
            if isempty(pool) || pool.NumWorkers~=o.ParallelWorkers
                if ~isempty(pool),delete(pool);end
                pool=parpool('Processes',double(o.ParallelWorkers));
            end
            active=parallel.FevalFuture.empty(0,1);
            activeRows=zeros(0,1);
            next=1;
            while numel(active)<pool.NumWorkers && next<=numel(pending)
                k=pending(next);
                active(end+1,1)=parfeval(pool,@one_restart,1, ...
                    taskMethod(k),taskRestart(k),cfg,xBaseline,masterSeed, ...
                    maxEvals,F.PatternSearch,signatureText,checkpointPaths{k},nativeDir);
                activeRows(end+1,1)=k;
                next=next+1;
            end
            while ~isempty(active)
                [idx,res]=fetchNext(active);
                k=activeRows(idx);
                active(idx)=[];
                activeRows(idx)=[];
                verify_run_result(res,taskMethod(k),taskRestart(k),Xfeasible(:,taskRestart(k)),cfg);
                results{k}=res;
                fprintf('Commit M %d/%d | %s restart %d | validated=%.9f\n', ...
                    nnz(~cellfun(@isempty,results)),height(taskPlan),char(res.method),res.restart,res.validatedAMI);
                if next<=numel(pending)
                    q=pending(next);
                    active(end+1,1)=parfeval(pool,@one_restart,1, ...
                        taskMethod(q),taskRestart(q),cfg,xBaseline,masterSeed, ...
                        maxEvals,F.PatternSearch,signatureText,checkpointPaths{q},nativeDir);
                    activeRows(end+1,1)=q;
                    next=next+1;
                end
            end
        else
            for ii=1:numel(pending)
                k=pending(ii);
                res=one_restart(taskMethod(k),taskRestart(k),cfg,xBaseline,masterSeed, ...
                    maxEvals,F.PatternSearch,signatureText,checkpointPaths{k},nativeDir);
                verify_run_result(res,taskMethod(k),taskRestart(k),Xfeasible(:,taskRestart(k)),cfg);
                results{k}=res;
            end
        end
    end

    assert(all(~cellfun(@isempty,results)),'CommitM:Incomplete','Not all 24 optimizer restarts completed.');
    Method=strings(24,1); Restart=zeros(24,1); Seed=zeros(24,1);
    ValidatedAMI=zeros(24,1); FastAMI=zeros(24,1);
    FunctionEvaluations=zeros(24,1); Iterations=zeros(24,1);
    Exitflag=nan(24,1); RuntimeSeconds=zeros(24,1);
    InitialLevels=strings(24,1); BestLevels=strings(24,1);
    for k=1:24
        r=results{k};
        Method(k)=r.method;Restart(k)=r.restart;Seed(k)=r.masterSeed;
        ValidatedAMI(k)=r.validatedAMI;FastAMI(k)=r.fastAMI;
        FunctionEvaluations(k)=r.functionEvaluations;
        Iterations(k)=r.iterations;Exitflag(k)=r.exitflag;
        RuntimeSeconds(k)=r.runtimeSeconds;
        InitialLevels(k)=string(mat2str(r.startFeasible(:).',12));
        BestLevels(k)=string(mat2str(r.xBest(:).',12));
    end
    restartTable=table(Method,Restart,Seed,ValidatedAMI,FastAMI, ...
        FunctionEvaluations,Iterations,Exitflag,RuntimeSeconds,InitialLevels,BestLevels);
    stats=commit_m_optimizer_statistics(restartTable);
    gpsBest=stats.BestValidatedAMI(stats.Method=="GPS");
    saBest=stats.BestValidatedAMI(stats.Method=="SA");
    four=table(names,[oldValues(1:2);gpsBest;saBest], ...
        'VariableNames',{'Method','ValidatedAMI'});
    assert(saBest<=log2(Msize)+1e-6 && gpsBest<=log2(Msize)+1e-6, ...
        'CommitM:AMIOutOfRange','AMI exceeds constellation entropy bound.');

    M=struct('version','CommitM-v1','dryRun',false, ...
        'runDirectory',runDir,'referenceCSVPath',refPath, ...
        'signature',signature,'taskPlan',taskPlan, ...
        'rawStarts',Xraw,'feasibleStarts',Xfeasible, ...
        'restartTable',restartTable,'optimizerStats',stats, ...
        'fourMethodTable',four,'results',{results}, ...
        'maxValidatedGapSAminusGPS',saBest-gpsBest, ...
        'aggregation','max over 12 per method; sample statistics over all 12');
    writetable(restartTable,fullfile(runDir,'CommitM_all_24_restarts.csv'));
    writetable(stats,fullfile(runDir,'CommitM_optimizer_stats_12_restarts.csv'));
    writetable(four,fullfile(runDir,'Fig7_CommitM_four_methods.csv'));
    M.figurePaths=plot_commit_m_fig7(four,runDir);
    save(resultFile,'M','-v7.3');
    fprintf('Commit M DONE: GPS max=%.9f SA max=%.9f gap=%+.9f\n',gpsBest,saBest,saBest-gpsBest);
    fprintf('Saved results, complete Figure 7 (PDF/PNG/FIG), and stats under %s\n',runDir);
end

function r=one_restart(method,restart,cfg,xBaseline,masterSeed,maxEvals,psPolicy,signatureText,checkpointPath,nativeDir)
    clock=tic;
    method=string(method);
    if method=="SA"
        s=sa_restart_task(cfg,xBaseline,masterSeed,restart,'[CommitM|SA]');
        r=struct('method',method,'restart',restart,'masterSeed',double(masterSeed), ...
            'startRaw',s.x0_raw(:),'startFeasible',s.x0_projected(:), ...
            'xBest',s.x_best(:),'fastAMI',s.bestMIFast, ...
            'validatedAMI',s.bestMIValidated, ...
            'functionEvaluations',cfg.SA.maxIter+1,'iterations',s.nIterations, ...
            'exitflag',NaN,'runtimeSeconds',toc(clock),'details',s);
    elseif method=="GPS"
        identity=struct('IdentityKey',sprintf('CommitM_GPS_r%02d',restart), ...
            'M',cfg.M,'SNRdB',cfg.SNR_dB,'SigmaR2',cfg.sigma_X_sq, ...
            'Replicate',1,'Restart',restart,'Seed',double(masterSeed));
        nativePath=fullfile(nativeDir,sprintf('CommitM_GPS_native_r%02d.mat',restart));
        g=i13d_patternsearch_restart_task(cfg,xBaseline,masterSeed,restart, ...
            maxEvals,psPolicy,identity,nativePath,signatureText);
        assert(g.success,'CommitM:GPSFailure','GPS restart %d failed: %s',restart,char(g.errorMessage));
        r=struct('method',method,'restart',restart,'masterSeed',double(masterSeed), ...
            'startRaw',g.x0Raw(:),'startFeasible',g.x0Feasible(:), ...
            'xBest',g.xBest(:),'fastAMI',g.bestMIFast, ...
            'validatedAMI',g.bestMIValidated, ...
            'functionEvaluations',g.functionEvaluations,'iterations',g.iterations, ...
            'exitflag',g.exitflag,'runtimeSeconds',toc(clock),'details',g);
    else
        error('CommitM:Method','Unsupported optimizer %s',method);
    end
    key=sprintf('%s_restart_%02d',char(method),restart);
    checkpoint=struct('signatureText',signatureText,'key',key,'result',r);
    temp=[tempname(fileparts(checkpointPath)) '.mat'];
    cleanup=onCleanup(@() delete_if_exists(temp)); %#ok<NASGU>
    save(temp,'checkpoint','-v7.3');
    [ok,msg]=movefile(temp,checkpointPath,'f');
    assert(ok,'CommitM:CheckpointSave','%s',msg);
end

function verify_run_result(r,method,restart,expectedX0,cfg)
    assert(string(r.method)==string(method) && r.restart==restart && ...
        numel(r.startFeasible)==numel(expectedX0) && ...
        norm(r.startFeasible(:)-expectedX0(:),inf)<1e-10 && ...
        isfinite(r.validatedAMI) && isfinite(r.fastAMI) && ...
        isfinite(r.functionEvaluations) && ...
        r.validatedAMI>=0 && r.validatedAMI<=log2(cfg.M)+1e-6 && ...
        r.functionEvaluations>=1 && ...
        abs(r.fastAMI-r.validatedAMI)<=1e-3, ...
        'CommitM:RunValidation','Optimizer result/start/AMI audit failed: %s restart %d',method,restart);
    AMI_functions.assert_constellation_feasible(r.xBest,cfg,1e-10);
end

function path=resolve_reference_csv(supplied)
    if strlength(string(supplied))>0
        path=char(supplied);
        assert(exist(path,'file')==2,'CommitM:MissingFrozenCSV','ReferenceCSVPath not found.');
        return;
    end
    root=fullfile(fso_result_utils.revision2_results_root(),'commit_k5p2_fig7');
    hits=dir(fullfile(root,'**','Fig7_K5p2_four_methods.csv'));
    assert(~isempty(hits),'CommitM:MissingFrozenCSV', ...
        'Supply ReferenceCSVPath pointing to the original K5.2 Fig7_K5p2_four_methods.csv.');
    [~,ix]=max([hits.datenum]);
    path=fullfile(hits(ix).folder,hits(ix).name);
end

function delete_if_exists(path)
    if exist(path,'file')==2
        try,delete(path);catch,end
    end
end

function tf=scalar_logical(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_int(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1&&v==fix(v);end
function tf=text_arg(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
