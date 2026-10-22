function Pfinal = run_commit_p_final_ieee_fig7(varargin)
%RUN_COMMIT_P_FINAL_IEEE_FIG7  Final IEEE Figure 7 run with new fixed seed.
% M8 / SNR20 / sigma_R^2=0.3; SA vs complete-poll GPS, MeshTolerance=1e-6.
% Exactly 12 paired unique starts for 24 runs. Separate baseline validators.
% Max over validated 12 values in four-bar no-errorbar grayscale figure.
% Save all-12 sample mean/variance/std, GPS termination and restart output.
% START/DONE from each worker with elapsed minutes and seconds.
%
%   addpath(genpath(fullfile(pwd,'src')));
%   test_commit_p_final_ieee_preflight();
%   plan=run_commit_p_final_ieee_fig7('DryRun',true);
%   Pfinal=run_commit_p_final_ieee_fig7();
%   Pfinal=run_commit_p_final_ieee_fig7('ResumeDirectory','/path/P_timestamp');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'M8ResultPath','',@is_text);
    addParameter(p,'M16ResultPath','',@is_text);
    addParameter(p,'OResultPath','',@is_text);
    addParameter(p,'ResumeDirectory','',@is_text);
    addParameter(p,'UseParallel',true,@is_bool);
    addParameter(p,'ParallelWorkers',8,@positive_int);
    addParameter(p,'DryRun',false,@is_bool);
    parse(p,varargin{:}); o=p.Results;
    P=i11p1_production_config();
    F=i13d_benchmark_config();
    F.PatternSearch.UseCompletePoll=true;
    assert(F.PatternSearch.MeshTolerance==1e-6 && ...
        F.PatternSearch.UseCompletePoll, ...
        'CommitP:GPSMesh','Expected exactly one mesh=1e-6 and CompletePoll=true.');
    m=8; snr=20; sigmaR2=0.3; nStarts=12; nTasks=2*nStarts;
    nSAIterations=4000; maxGPSEvals=nSAIterations+1;
    [loghDt,dtMeta]=select_logh_dt(m,snr,sigmaR2, ...
        'Policy',P.LoghDtPolicy,'BaseDt',P.loghDt, ...
        'RefinedDt',P.loghDtRefined);
    assert(loghDt==0.01 && F.MaxEvalsPerStart(F.MVec==m)==maxGPSEvals && ...
        P.P_avg==1 && P.minGap==0.01 && P.pinZero && ...
        P.T0==0.4 && P.Tf==1e-3 && P.BaseStd0==0.15 && ...
        P.ItersPerTemp==50 && strcmp(P.FastFadingMethod,'logh') && ...
        strcmp(P.YGridMode,'candidate-adaptive') && P.YGridScale==1.1 && ...
        P.loghSpanSigma==8 && P.loghYBlockSize==512 && ...
        F.PatternSearch.UseCompletePoll && ~F.PatternSearch.UseParallel, ...
        'CommitP:FrozenConfig','The scientific conditions differ from Commit M.');
    cfg=build_fso_config(m,P.P_avg,snr,sigmaR2, ...
        'FastFadingMethod',P.FastFadingMethod,'loghDt',loghDt, ...
        'loghSpanSigma',P.loghSpanSigma, ...
        'loghYBlockSize',P.loghYBlockSize, ...
        'YGridMode',P.YGridMode,'YGridScale',P.YGridScale, ...
        'minGap',P.minGap,'pinZero',P.pinZero, ...
        'saMaxIter',nSAIterations,'saNStarts',nStarts, ...
        'seedInit',0,'logEvery',0,'historyEvery',P.historyEvery);
    cfg.SA.T0=P.T0;cfg.SA.Tf=P.Tf;
    cfg.SA.baseStd0=P.BaseStd0;cfg.SA.itersPerTemp=P.ItersPerTemp;
    cfg.SA.nBlocks=ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate=exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);

    uniform=define_constellation(m,P.P_avg,true,"mean");
    % Commit P pre-declares a new deterministic seed, distinct from M.
    seedOffset=104729;
    masterSeed=fso_result_utils.case_seed(P.baseSeed+seedOffset,m,snr,sigmaR2);
    commitMSeed=fso_result_utils.case_seed(P.baseSeed,m,snr,sigmaR2);
    assert(masterSeed~=commitMSeed,'CommitP:Seed','New seed must differ from M.');
    Xraw=i13_shared_start_points(cfg,uniform,masterSeed,nStarts);
    Xfeasible=zeros(m,nStarts);
    for k=1:nStarts
        Xfeasible(:,k)=AMI_functions.project_constellation_1D(Xraw(:,k),cfg);
        AMI_functions.assert_constellation_feasible(Xfeasible(:,k),cfg,1e-10);
    end
    assert(norm(Xfeasible(:,1)-uniform(:),inf)<1e-12, ...
        'CommitP:FirstStart','First restart must be uniform PAM.');
    for k=1:nStarts
        for j=k+1:nStarts
            assert(norm(Xfeasible(:,k)-Xfeasible(:,j),inf)>1e-10, ...
                'CommitP:DuplicateStart','Restart starting points must be distinct.');
        end
    end
    taskMethod=[repmat("SA",nStarts,1);repmat("GPS",nStarts,1)];
    taskRestart=repmat((1:nStarts)',2,1);
    taskPlan=table(taskMethod,taskRestart,'VariableNames',{'Method','Restart'});
    frozen=struct('version','CommitP-v1','M',m,'SNRdB',snr, ...
        'SigmaR2',sigmaR2,'Pavg',P.P_avg,'minGap',P.minGap, ...
        'nStarts',nStarts,'masterSeed',masterSeed, ...
        'seedOffset',seedOffset,'commitMSeed',commitMSeed, ...
        'SAIterations',nSAIterations,'GPSMaxEvals',maxGPSEvals, ...
        'saT0',P.T0,'saTf',P.Tf,'saBaseStd0',P.BaseStd0, ...
        'saItersPerTemp',P.ItersPerTemp,'dtMeta',dtMeta, ...
        'loghDt',loghDt,'loghSpanSigma',P.loghSpanSigma, ...
        'YGridMode',P.YGridMode,'YGridScale',P.YGridScale, ...
        'loghYBlockSize',P.loghYBlockSize, ...
        'GPSOptions',F.PatternSearch,'ZhouExtraBits',2, ...
        'startPolicy','uniform first, 11 Threefry gaussian 0.2, shared SA/GPS');
    signature=char(jsonencode(frozen));

    if o.DryRun
        Pfinal=struct('dryRun',true,'config',frozen,'taskPlan',taskPlan, ...
            'rawStarts',Xraw,'feasibleStarts',Xfeasible);
        fprintf('Commit P DRY RUN: M8 SNR20, 12 distinct paired starts, 24 tasks, GPS mesh=1e-6.\n');
        return;
    end

    pathM=find_m8(o.M8ResultPath);
    q=load(pathM,'M');
    assert(isfield(q,'M') && isfield(q.M,'signature') && ...
        isfield(q.M,'fourMethodTable') && isfield(q.M,'optimizerStats') && ...
        q.M.signature.M==8 && q.M.signature.SNRdB==20 && ...
        abs(q.M.signature.SigmaR2-sigmaR2)<1e-12 && ...
        q.M.signature.nStarts==12 && q.M.signature.PatternSearch.UseCompletePoll, ...
        'CommitP:HistoryM','Expected finished Commit M at M8 SNR20.');
    oldM=q.M;
    pathN=find_m16(o.M16ResultPath);
    q=load(pathN,'N');
    assert(isfield(q,'N') && isfield(q.N,'config') && ...
        isfield(q.N,'fourMethodTable') && isfield(q.N,'optimizerStats') && ...
        q.N.config.M==16 && q.N.config.SNRdB==20 && ...
        abs(q.N.config.SigmaR2-sigmaR2)<1e-12 && ...
        q.N.config.nStarts==12 && q.N.config.GPSOptions.UseCompletePoll, ...
        'CommitP:HistoryN','Expected finished Commit N at M16 SNR20.');
    oldN=q.N;
    pathO=find_o(o.OResultPath);
    q=load(pathO,'O');
    assert(isfield(q,'O') && isfield(q.O,'config') && ...
        isfield(q.O,'fourMethodTable') && isfield(q.O,'optimizerStats') && ...
        q.O.config.M==8 && q.O.config.SNRdB==30 && ...
        abs(q.O.config.SigmaR2-sigmaR2)<1e-12 && ...
        q.O.config.nStarts==12 && q.O.config.GPSOptions.UseCompletePoll, ...
        'CommitP:HistoryO','Expected finished Commit O at M8 SNR30.');
    oldO=q.O;

    % Independently recompute both baselines at the final M8 SNR20 point.
    fprintf('[CommitP|M8|SNR20|Uniform PAM] START\n');
    t=tic; uniformMI=cfg.AMI_Validator(uniform(:));
    fprintf('[CommitP|M8|SNR20|Uniform PAM] DONE | AMI=%.9f | elapsed=%s\n',uniformMI,mmss(toc(t)));
    fprintf('[CommitP|M8|SNR20|Zhou GS] START\n');
    t=tic;
    zhou=zhou2025_optical_gs_constellation(m,P.P_avg, ...
        'ExtraBits',2,'Regularize',true);
    AMI_functions.assert_constellation_feasible(zhou.x,cfg,1e-10);
    zhouMI=cfg.AMI_Validator(zhou.x);
    fprintf('[CommitP|M8|SNR20|Zhou GS] DONE | AMI=%.9f | elapsed=%s\n',zhouMI,mmss(toc(t)));
    assert(all(isfinite([uniformMI,zhouMI])) && ...
        all([uniformMI,zhouMI]>=0) && all([uniformMI,zhouMI]<=log2(m)+1e-6), ...
        'CommitP:BadBaselines','Invalid baseline AMI.');

    outputRoot=fullfile(fso_result_utils.revision2_results_root(),'commit_p_final_ieee_fig7');
    if strlength(string(o.ResumeDirectory))>0
        outDir=char(o.ResumeDirectory);
        assert(isfolder(outDir),'CommitP:BadResume','ResumeDirectory does not exist.');
    else
        stamp=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
        outDir=fullfile(outputRoot,['P_' stamp]);
        if ~isfolder(outDir),mkdir(outDir);end
    end
    cpDir=fullfile(outDir,'checkpoints');
    gpsDir=fullfile(outDir,'native_gps');
    if ~isfolder(cpDir),mkdir(cpDir);end
    if ~isfolder(gpsDir),mkdir(gpsDir);end

    results=cell(nTasks,1); cp=cell(nTasks,1);
    for k=1:nTasks
        key=sprintf('%s_restart_%02d',char(taskMethod(k)),taskRestart(k));
        cp{k}=fullfile(cpDir,[key '.mat']);
        if isfile(cp{k})
            L=load(cp{k},'checkpoint');
            assert(isfield(L,'checkpoint') && ...
                strcmp(L.checkpoint.signature,signature) && ...
                strcmp(L.checkpoint.key,key), ...
                'CommitP:Checkpoint','Checkpoint science/identity mismatch: %s',cp{k});
            results{k}=L.checkpoint.result;
            audit_one(results{k},taskMethod(k),taskRestart(k), ...
                Xfeasible(:,taskRestart(k)),cfg,maxGPSEvals);
            fprintf('[CommitP|M8|SNR20|%s|restart %02d] RESUMED | validated=%.9f | elapsed=%s\n', ...
                char(taskMethod(k)),taskRestart(k),results{k}.validatedAMI,mmss(results{k}.runtimeSeconds));
        end
    end

    pending=find(cellfun(@isempty,results));
    fprintf('Commit P: %d/24 completed, %d pending | %s\n',nTasks-numel(pending),numel(pending),outDir);
    if ~isempty(pending)
        progress=parallel.pool.DataQueue;
        afterEach(progress,@(message) fprintf('%s\n',message));
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
                row=pending(next);
                active(end+1,1)=parfeval(pool,@run_one,1, ...
                    taskMethod(row),taskRestart(row),cfg,uniform,masterSeed, ...
                    maxGPSEvals,F.PatternSearch,signature,cp{row},gpsDir,progress);
                activeRows(end+1,1)=row;
                next=next+1;
            end
            try
                while ~isempty(active)
                    [finished,r]=fetchNext(active);
                    row=activeRows(finished);
                    active(finished)=[];
                    activeRows(finished)=[];
                    audit_one(r,taskMethod(row),taskRestart(row), ...
                        Xfeasible(:,taskRestart(row)),cfg,maxGPSEvals);
                    results{row}=r;
                    fprintf('[CommitP|M8|SNR20] Progress: %d/%d restarts completed\n', ...
                        nnz(~cellfun(@isempty,results)),nTasks);
                    if next<=numel(pending)
                        newer=pending(next);
                        active(end+1,1)=parfeval(pool,@run_one,1, ...
                            taskMethod(newer),taskRestart(newer),cfg,uniform,masterSeed, ...
                            maxGPSEvals,F.PatternSearch,signature,cp{newer},gpsDir,progress);
                        activeRows(end+1,1)=newer;
                        next=next+1;
                    end
                end
            catch ME
                if ~isempty(active),cancel(active);end
                rethrow(ME);
            end
        else
            for ii=1:numel(pending)
                row=pending(ii);
                results{row}=run_one(taskMethod(row),taskRestart(row), ...
                    cfg,uniform,masterSeed,maxGPSEvals,F.PatternSearch, ...
                    signature,cp{row},gpsDir,progress);
                audit_one(results{row},taskMethod(row),taskRestart(row), ...
                    Xfeasible(:,taskRestart(row)),cfg,maxGPSEvals);
            end
        end
    end

    assert(all(~cellfun(@isempty,results)),'CommitP:Incomplete','24 validated tasks are required.');
    Method=strings(nTasks,1);Restart=zeros(nTasks,1);
    ValidatedAMI=zeros(nTasks,1);FastAMI=zeros(nTasks,1);
    FunctionEvaluations=zeros(nTasks,1);Iterations=zeros(nTasks,1);
    Exitflag=nan(nTasks,1);RuntimeSeconds=zeros(nTasks,1);
    FinalMeshSize=nan(nTasks,1);MeshTolerance=nan(nTasks,1);
    StartLevels=strings(nTasks,1);BestLevels=strings(nTasks,1);
    for k=1:nTasks
        r=results{k};Method(k)=r.method;Restart(k)=r.restart;
        ValidatedAMI(k)=r.validatedAMI;FastAMI(k)=r.fastAMI;
        FunctionEvaluations(k)=r.functionEvaluations;
        Iterations(k)=r.iterations;Exitflag(k)=r.exitflag;
        RuntimeSeconds(k)=r.runtimeSeconds;
        if r.method=="GPS"
            FinalMeshSize(k)=r.details.finalMeshSize;
            MeshTolerance(k)=F.PatternSearch.MeshTolerance;
        end
        StartLevels(k)=string(mat2str(r.startFeasible(:).',12));
        BestLevels(k)=string(mat2str(r.xBest(:).',12));
    end
    restartTable=table(Method,Restart,ValidatedAMI,FastAMI, ...
        FunctionEvaluations,Iterations,Exitflag,RuntimeSeconds, ...
        FinalMeshSize,MeshTolerance,StartLevels,BestLevels);
    stats=commit_m_optimizer_statistics(restartTable);
    gpsMax=stats.BestValidatedAMI(stats.Method=="GPS");
    saMax=stats.BestValidatedAMI(stats.Method=="SA");

    methods=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    four=table(methods,[uniformMI;zhouMI;gpsMax;saMax], ...
        'VariableNames',{'Method','ValidatedAMI'});
    methodsExpected=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    tableM=oldM.fourMethodTable;tableN=oldN.fourMethodTable;
    tableO=oldO.fourMethodTable;
    assert(isequal(string(four.Method(:)),methodsExpected) && ...
        isequal(string(tableM.Method(:)),methodsExpected) && ...
        isequal(string(tableN.Method(:)),methodsExpected) && ...
        isequal(string(tableO.Method(:)),methodsExpected), ...
        'CommitP:MethodOrder','Previous method ordering mismatch.');
    im=double(tableM.ValidatedAMI(:));in=double(tableN.ValidatedAMI(:));
    io=double(tableO.ValidatedAMI(:));ip=double(four.ValidatedAMI(:));
    historical=table(methodsExpected,im,in,io,ip, ...
        im/log2(8),in/log2(16),io/log2(8),ip/log2(8), ...
        'VariableNames',{'Method','CommitM_M8_SNR20','CommitN_M16_SNR20', ...
        'CommitO_M8_SNR30','CommitP_M8_SNR20', ...
        'NormalizedM','NormalizedN','NormalizedO','NormalizedP'});
    a=oldM.optimizerStats;b=oldN.optimizerStats;c=oldO.optimizerStats;
    c=c(ismember(string(c.Method),["GPS","SA"]),:);
    mNames=["GPS";"SA"];
    assert(isequal(string(a.Method(:)),mNames) && ...
        isequal(string(b.Method(:)),mNames) && ...
        isequal(string(c.Method(:)),mNames) && ...
        isequal(string(stats.Method(:)),mNames), ...
        'CommitP:OptimizerStatsOrder','Previous optimizer order mismatch.');
    historicalStats=table(mNames, ...
        a.BestValidatedAMI,b.BestValidatedAMI,c.BestValidatedAMI,stats.BestValidatedAMI, ...
        a.MeanValidatedAMI,b.MeanValidatedAMI,c.MeanValidatedAMI,stats.MeanValidatedAMI, ...
        a.StdValidatedAMI,b.StdValidatedAMI,c.StdValidatedAMI,stats.StdValidatedAMI, ...
        'VariableNames',{'Method','MMax','NMax','OMax','PMax', ...
        'MMean','NMean','OMean','PMean','MStd','NStd','OStd','PStd'});
    repeatability=table(methodsExpected,im,ip,ip-im, ...
        'VariableNames',{'Method','CommitMValidatedAMI', ...
        'CommitPValidatedAMI','CommitPMinusM'});
    Pfinal=struct('version','CommitP-v1','dryRun',false,'config',frozen, ...
        'outputDirectory',outDir,'sourceCommitM',pathM,'sourceCommitN',pathN, ...
        'sourceCommitO',pathO,'taskPlan',taskPlan, ...
        'rawStarts',Xraw,'feasibleStarts',Xfeasible, ...
        'uniformConstellation',uniform(:),'zhouGS',zhou, ...
        'restartTable',restartTable,'optimizerStats',stats, ...
        'fourMethodTable',four,'historicalComparison',historical, ...
        'historicalOptimizerStats',historicalStats, ...
        'repeatabilityCommitM',repeatability,'results',{results});
    writetable(restartTable,fullfile(outDir,'CommitP_all_24_restarts.csv'));
    writetable(stats,fullfile(outDir,'CommitP_optimizer_stats.csv'));
    writetable(four,fullfile(outDir,'Fig7_CommitP_four_methods.csv'));
    writetable(historical,fullfile(outDir,'CommitP_comparison_M_N_O.csv'));
    writetable(historicalStats,fullfile(outDir,'CommitP_optimizer_comparison_M_N_O.csv'));
    writetable(repeatability,fullfile(outDir,'CommitP_repeatability_vs_M.csv'));
    Pfinal.figurePaths=plot_commit_p_final_ieee_fig7(four,outDir);
    save(fullfile(outDir,'CommitP_results.mat'),'Pfinal','-v7.3');
    fprintf('Commit P DONE | GPS=%.9f SA=%.9f gap=%+.9f\n', ...
        gpsMax,saMax,saMax-gpsMax);
    disp(historical);disp(repeatability);
    fprintf('All outputs: %s\n',outDir);
end

function r=run_one(method,k,cfg,uniform,seed,maxEvals,policy,signature,cpPath,gpsDir,progress)
    method=string(method);
    timer=tic;
    send(progress,sprintf('[CommitP|M8|SNR20|%s|restart %02d] START',char(method),k));
    try
        if method=="SA"
            a=sa_restart_task(cfg,uniform,seed,k,'[CommitP|SA]');
            r=struct('method',method,'restart',k,'startRaw',a.x0_raw(:), ...
                'startFeasible',a.x0_projected(:),'xBest',a.x_best(:), ...
                'fastAMI',a.bestMIFast,'validatedAMI',a.bestMIValidated, ...
                'functionEvaluations',cfg.SA.maxIter+1, ...
                'iterations',a.nIterations,'exitflag',NaN,'details',a);
        elseif method=="GPS"
            assert(policy.MeshTolerance==1e-6 && policy.UseCompletePoll, ...
                'CommitP:GPSPolicy','GPS policy must remain complete and mesh=1e-6.');
            ident=struct('IdentityKey',sprintf('CommitP_%s_r%02d',char(method),k), ...
                'M',cfg.M,'SNRdB',cfg.SNR_dB,'SigmaR2',cfg.sigma_X_sq, ...
                'Replicate',1,'Restart',k,'Seed',double(seed));
            native=fullfile(gpsDir,sprintf('%s_r%02d.mat',char(method),k));
            g=i13d_patternsearch_restart_task(cfg,uniform,seed,k,maxEvals, ...
                policy,ident,native,signature);
            assert(g.success,'CommitP:GPSFailure','%s restart %d: %s', ...
                char(method),k,char(g.errorMessage));
            r=struct('method',method,'restart',k,'startRaw',g.x0Raw(:), ...
                'startFeasible',g.x0Feasible(:),'xBest',g.xBest(:), ...
                'fastAMI',g.bestMIFast,'validatedAMI',g.bestMIValidated, ...
                'functionEvaluations',g.functionEvaluations, ...
                'iterations',g.iterations,'exitflag',g.exitflag,'details',g);
        else
            error('CommitP:Method','Unknown method %s',char(method));
        end
        r.runtimeSeconds=toc(timer);
        key=sprintf('%s_restart_%02d',char(method),k);
        checkpoint=struct('signature',signature,'key',key,'result',r);
        temp=[tempname(fileparts(cpPath)) '.mat'];
        cleanup=onCleanup(@()delete_tmp(temp)); %#ok<NASGU>
        save(temp,'checkpoint','-v7.3');
        [ok,msg]=movefile(temp,cpPath,'f');
        assert(ok,'CommitP:CheckpointSave','%s',msg);
        send(progress,sprintf('[CommitP|M8|SNR20|%s|restart %02d] DONE | elapsed=%s | validated=%.9f | evals=%d', ...
            char(method),k,mmss(r.runtimeSeconds),r.validatedAMI,r.functionEvaluations));
    catch ME
        send(progress,sprintf('[CommitP|M8|SNR20|%s|restart %02d] FAILED | elapsed=%s | %s', ...
            char(method),k,mmss(toc(timer)),ME.message));
        rethrow(ME);
    end
end

function audit_one(r,method,k,x0,cfg,budget)
    assert(string(r.method)==string(method) && r.restart==k && ...
        norm(r.startFeasible(:)-x0(:),inf)<1e-10 && ...
        isfinite(r.validatedAMI) && r.validatedAMI>=0 && ...
        r.validatedAMI<=log2(cfg.M)+1e-6 && ...
        abs(r.fastAMI-r.validatedAMI)<1e-3 && ...
        isfinite(r.functionEvaluations) && ...
        r.functionEvaluations>=1 && r.functionEvaluations<=budget, ...
        'CommitP:ResultAudit','Run result invalid for %s restart %d.',method,k);
    AMI_functions.assert_constellation_feasible(r.xBest,cfg,1e-10);
end

function path=find_m8(provided)
    if strlength(string(provided))>0
        path=char(provided);
        assert(isfile(path),'CommitP:M8Path','Provided M8ResultPath not found.');
        return;
    end
    root=fullfile(fso_result_utils.revision2_results_root(),'commit_m_final_fig7');
    hits=dir(fullfile(root,'**','CommitM_results.mat'));
    assert(~isempty(hits),'CommitP:MissingM8', ...
        'Please pass M8ResultPath pointing to your finished CommitM_results.mat.');
    [~,ii]=max([hits.datenum]);
    path=fullfile(hits(ii).folder,hits(ii).name);
end

function path=find_m16(provided)
    if strlength(string(provided))>0
        path=char(provided);
        assert(isfile(path),'CommitP:MissingM16','Provided M16ResultPath not found.');
        return;
    end
    root=fullfile(fso_result_utils.revision2_results_root(),'commit_n_m16_fig7');
    hits=dir(fullfile(root,'**','CommitN_M16_results.mat'));
    assert(~isempty(hits),'CommitP:MissingM16', ...
        'Pass M16ResultPath for the completed Commit N result.');
    [~,ii]=max([hits.datenum]);
    path=fullfile(hits(ii).folder,hits(ii).name);
end

function path=find_o(provided)
    if strlength(string(provided))>0
        path=char(provided);
        assert(isfile(path),'CommitP:MissingO','OResultPath not found.');
        return;
    end
    root=fullfile(fso_result_utils.revision2_results_root(),'commit_o_m8_snr30');
    hits=dir(fullfile(root,'**','CommitO_results.mat'));
    assert(~isempty(hits),'CommitP:MissingO','Pass OResultPath for completed Commit O.');
    [~,ii]=max([hits.datenum]);
    path=fullfile(hits(ii).folder,hits(ii).name);
end

function str=mmss(seconds)
    seconds=max(0,double(seconds));
    str=sprintf('%02dm %05.2fs',floor(seconds/60),mod(seconds,60));
end

function delete_tmp(file)
    if isfile(file)
        try,delete(file);catch,end
    end
end
function tf=is_text(x),tf=ischar(x)||(isstring(x)&&isscalar(x));end
function tf=is_bool(x),tf=islogical(x)&&isscalar(x);end
function tf=positive_int(x),tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x);end
