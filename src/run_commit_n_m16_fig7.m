function N = run_commit_n_m16_fig7(varargin)
%RUN_COMMIT_N_M16_FIG7  M16 paired SA/GPS, reference baselines, and M8 comparison.
% 12 unique matched raw and projected starting points, same frozen Commit M
% settings except M=16. GPS uses complete polling. Each method reports the
% maximum independently validated AMI over 12 restarts; mean, sample
% variance and sample std use all 12. Figure has NO error bars.
%
% From the repository root:
%   addpath(genpath(fullfile(pwd,'src')));
%   test_commit_n_m16_preflight();
%   plan=run_commit_n_m16_fig7('DryRun',true);
%   N=run_commit_n_m16_fig7();
%
% To override M8 source:
%   N=run_commit_n_m16_fig7('M8ResultPath','/path/CommitM_results.mat');
% To resume:
%   N=run_commit_n_m16_fig7('ResumeDirectory','/path/N_timestamp');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'M8ResultPath','',@is_text);
    addParameter(p,'ResumeDirectory','',@is_text);
    addParameter(p,'UseParallel',true,@is_bool);
    addParameter(p,'ParallelWorkers',8,@positive_int);
    addParameter(p,'DryRun',false,@is_bool);
    parse(p,varargin{:}); o=p.Results;
    P=i11p1_production_config();
    F=i13d_benchmark_config();
    F.PatternSearch.UseCompletePoll=true;
    m=16; snr=20; sigmaR2=0.3; nStarts=12; nTasks=2*nStarts;
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
        'CommitN:FrozenConfig','The scientific conditions differ from Commit M.');
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
    masterSeed=fso_result_utils.case_seed(P.baseSeed,m,snr,sigmaR2);
    Xraw=i13_shared_start_points(cfg,uniform,masterSeed,nStarts);
    Xfeasible=zeros(m,nStarts);
    for k=1:nStarts
        Xfeasible(:,k)=AMI_functions.project_constellation_1D(Xraw(:,k),cfg);
        AMI_functions.assert_constellation_feasible(Xfeasible(:,k),cfg,1e-10);
    end
    assert(norm(Xfeasible(:,1)-uniform(:),inf)<1e-12, ...
        'CommitN:FirstStart','First restart must be uniform PAM.');
    for k=1:nStarts
        for j=k+1:nStarts
            assert(norm(Xfeasible(:,k)-Xfeasible(:,j),inf)>1e-10, ...
                'CommitN:DuplicateStart','Restart starting points must be distinct.');
        end
    end
    taskMethod=[repmat("SA",nStarts,1);repmat("GPS",nStarts,1)];
    taskRestart=[(1:nStarts)';(1:nStarts)'];
    taskPlan=table(taskMethod,taskRestart,'VariableNames',{'Method','Restart'});
    frozen=struct('version','CommitN-v1','M',m,'SNRdB',snr, ...
        'SigmaR2',sigmaR2,'Pavg',P.P_avg,'minGap',P.minGap, ...
        'nStarts',nStarts,'masterSeed',masterSeed, ...
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
        N=struct('dryRun',true,'config',frozen,'taskPlan',taskPlan, ...
            'rawStarts',Xraw,'feasibleStarts',Xfeasible);
        fprintf('Commit N DRY RUN: M=16, 12 unique paired starts, 24 tasks, complete GPS poll.\n');
        return;
    end

    % Load and audit the previous M=8 result BEFORE running expensive work.
    m8Path=find_m8(o.M8ResultPath);
    C=load(m8Path,'M');
    assert(isfield(C,'M') && isfield(C.M,'signature') && ...
        isfield(C.M,'optimizerStats') && isfield(C.M,'fourMethodTable') && ...
        C.M.signature.M==8 && C.M.signature.SNRdB==snr && ...
        abs(C.M.signature.SigmaR2-sigmaR2)<1e-12 && ...
        C.M.signature.nStarts==12 && ...
        C.M.signature.PatternSearch.UseCompletePoll, ...
        'CommitN:InvalidM8','Expected complete Commit M result for M8/20dB/0.3.');

    % Independently evaluate BOTH baselines for the new M=16 constellation.
    fprintf('[CommitN|M16|Uniform PAM] START\n');
    t=tic; uniformMI=cfg.AMI_Validator(uniform(:));
    fprintf('[CommitN|M16|Uniform PAM] DONE | AMI=%.9f | elapsed=%.2fs\n',uniformMI,toc(t));
    fprintf('[CommitN|M16|Zhou GS] START\n');
    t=tic;
    zhou=zhou2025_optical_gs_constellation(m,P.P_avg, ...
        'ExtraBits',2,'Regularize',true);
    AMI_functions.assert_constellation_feasible(zhou.x,cfg,1e-10);
    zhouMI=cfg.AMI_Validator(zhou.x);
    fprintf('[CommitN|M16|Zhou GS] DONE | AMI=%.9f | elapsed=%.2fs\n',zhouMI,toc(t));
    assert(all(isfinite([uniformMI,zhouMI])) && ...
        all([uniformMI,zhouMI]>=0) && all([uniformMI,zhouMI]<=log2(m)+1e-6), ...
        'CommitN:BadBaselines','Invalid baseline AMI.');

    outputRoot=fullfile(fso_result_utils.revision2_results_root(),'commit_n_m16_fig7');
    if strlength(string(o.ResumeDirectory))>0
        outDir=char(o.ResumeDirectory);
        assert(isfolder(outDir),'CommitN:BadResume','ResumeDirectory does not exist.');
    else
        stamp=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
        outDir=fullfile(outputRoot,['N_' stamp]);
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
                'CommitN:Checkpoint','Checkpoint science/identity mismatch: %s',cp{k});
            results{k}=L.checkpoint.result;
            audit_one(results{k},taskMethod(k),taskRestart(k), ...
                Xfeasible(:,taskRestart(k)),cfg,maxGPSEvals);
            fprintf('[CommitN|M16|%s|restart %02d] RESUMED | validated=%.9f | elapsed=%.2fs\n', ...
                char(taskMethod(k)),taskRestart(k),results{k}.validatedAMI,results{k}.runtimeSeconds);
        end
    end

    pending=find(cellfun(@isempty,results));
    fprintf('Commit N: %d/24 completed, %d pending | %s\n',nTasks-numel(pending),numel(pending),outDir);
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
                    fprintf('[CommitN|M16] Progress: %d/%d restarts completed\n', ...
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

    assert(all(~cellfun(@isempty,results)),'CommitN:Incomplete','24 validated tasks are required.');
    Method=strings(nTasks,1);Restart=zeros(nTasks,1);
    ValidatedAMI=zeros(nTasks,1);FastAMI=zeros(nTasks,1);
    FunctionEvaluations=zeros(nTasks,1);Iterations=zeros(nTasks,1);
    Exitflag=nan(nTasks,1);RuntimeSeconds=zeros(nTasks,1);
    StartLevels=strings(nTasks,1);BestLevels=strings(nTasks,1);
    for k=1:nTasks
        r=results{k};Method(k)=r.method;Restart(k)=r.restart;
        ValidatedAMI(k)=r.validatedAMI;FastAMI(k)=r.fastAMI;
        FunctionEvaluations(k)=r.functionEvaluations;
        Iterations(k)=r.iterations;Exitflag(k)=r.exitflag;
        RuntimeSeconds(k)=r.runtimeSeconds;
        StartLevels(k)=string(mat2str(r.startFeasible(:).',12));
        BestLevels(k)=string(mat2str(r.xBest(:).',12));
    end
    restartTable=table(Method,Restart,ValidatedAMI,FastAMI, ...
        FunctionEvaluations,Iterations,Exitflag,RuntimeSeconds,StartLevels,BestLevels);
    stats=commit_m_optimizer_statistics(restartTable);
    gpsMax=stats.BestValidatedAMI(stats.Method=="GPS");
    saMax=stats.BestValidatedAMI(stats.Method=="SA");

    methods=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    four=table(methods,[uniformMI;zhouMI;gpsMax;saMax], ...
        'VariableNames',{'Method','ValidatedAMI'});
    old=C.M.fourMethodTable;
    assert(istable(old)&&height(old)==4 && ...
        isequal(string(old.Method(:)),methods), ...
        'CommitN:M8Schema','M8 four-method data not compatible.');
    I8=double(old.ValidatedAMI(:)); I16=double(four.ValidatedAMI(:));
    comparison=table(methods,I8,I16,I16-I8,I8/log2(8),I16/log2(16), ...
        'VariableNames',{'Method','M8ValidatedAMI','M16ValidatedAMI', ...
        'DeltaBitsPerSymbol','M8FractionOfLog2M','M16FractionOfLog2M'});
    S8=C.M.optimizerStats;
    assert(isequal(string(S8.Method(:)),string(stats.Method(:))), ...
        'CommitN:M8StatsSchema','M8 optimizer statistics missing or reordered.');
    comparisonStats=table(string(stats.Method), ...
        double(S8.BestValidatedAMI),double(stats.BestValidatedAMI), ...
        double(S8.MeanValidatedAMI),double(stats.MeanValidatedAMI), ...
        double(S8.StdValidatedAMI),double(stats.StdValidatedAMI), ...
        'VariableNames',{'Method','M8Max','M16Max', ...
        'M8Mean','M16Mean','M8Std12','M16Std12'});
    N=struct('version','CommitN-v1','dryRun',false,'config',frozen, ...
        'outputDirectory',outDir,'m8ResultPath',m8Path, ...
        'taskPlan',taskPlan,'rawStarts',Xraw,'feasibleStarts',Xfeasible, ...
        'uniformConstellation',uniform(:),'zhouGS',zhou, ...
        'restartTable',restartTable,'optimizerStats',stats, ...
        'fourMethodTable',four,'m8VsM16',comparison, ...
        'm8VsM16OptimizerStats',comparisonStats,'results',{results});
    writetable(restartTable,fullfile(outDir,'CommitN_M16_all_24_restarts.csv'));
    writetable(stats,fullfile(outDir,'CommitN_M16_optimizer_stats_12_restarts.csv'));
    writetable(four,fullfile(outDir,'Fig7_CommitN_M16_four_methods.csv'));
    writetable(comparison,fullfile(outDir,'CommitN_M8_vs_M16_four_methods.csv'));
    writetable(comparisonStats,fullfile(outDir,'CommitN_M8_vs_M16_optimizer_stats.csv'));
    N.figurePaths=plot_commit_n_m16_fig7(four,outDir);
    save(fullfile(outDir,'CommitN_M16_results.mat'),'N','-v7.3');
    fprintf('Commit N DONE | M=16 GPS=%.9f SA=%.9f | SA-GPS=%+.9f\n',gpsMax,saMax,saMax-gpsMax);
    disp(comparison);
    fprintf('Outputs saved under %s\n',outDir);
end

function r=run_one(method,k,cfg,uniform,seed,maxEvals,policy,signature,cpPath,gpsDir,progress)
    method=string(method);
    timer=tic;
    send(progress,sprintf('[CommitN|M16|%s|restart %02d] START',char(method),k));
    try
        if method=="SA"
            s=sa_restart_task(cfg,uniform,seed,k,'[CommitN|SA]');
            r=struct('method',method,'restart',k,'startRaw',s.x0_raw(:), ...
                'startFeasible',s.x0_projected(:),'xBest',s.x_best(:), ...
                'fastAMI',s.bestMIFast,'validatedAMI',s.bestMIValidated, ...
                'functionEvaluations',cfg.SA.maxIter+1, ...
                'iterations',s.nIterations,'exitflag',NaN,'details',s);
        else
            ident=struct('IdentityKey',sprintf('CommitN_M16_GPS_r%02d',k), ...
                'M',cfg.M,'SNRdB',cfg.SNR_dB,'SigmaR2',cfg.sigma_X_sq, ...
                'Replicate',1,'Restart',k,'Seed',double(seed));
            gpsCheckpoint=fullfile(gpsDir,sprintf('GPS_r%02d.mat',k));
            g=i13d_patternsearch_restart_task(cfg,uniform,seed,k,maxEvals, ...
                policy,ident,gpsCheckpoint,signature);
            assert(g.success,'CommitN:GPSFailed','GPS %d failed: %s',k,char(g.errorMessage));
            r=struct('method',method,'restart',k,'startRaw',g.x0Raw(:), ...
                'startFeasible',g.x0Feasible(:),'xBest',g.xBest(:), ...
                'fastAMI',g.bestMIFast,'validatedAMI',g.bestMIValidated, ...
                'functionEvaluations',g.functionEvaluations, ...
                'iterations',g.iterations,'exitflag',g.exitflag,'details',g);
        end
        r.runtimeSeconds=toc(timer);
        key=sprintf('%s_restart_%02d',char(method),k);
        checkpoint=struct('signature',signature,'key',key,'result',r);
        tmp=[tempname(fileparts(cpPath)) '.mat'];
        cleanup=onCleanup(@()delete_tmp(tmp)); %#ok<NASGU>
        save(tmp,'checkpoint','-v7.3');
        [ok,msg]=movefile(tmp,cpPath,'f');
        assert(ok,'CommitN:CheckpointSave','%s',msg);
        send(progress,sprintf('[CommitN|M16|%s|restart %02d] DONE | elapsed=%.2fs | validated=%.9f | evals=%d', ...
            char(method),k,r.runtimeSeconds,r.validatedAMI,r.functionEvaluations));
    catch ME
        send(progress,sprintf('[CommitN|M16|%s|restart %02d] FAILED | elapsed=%.2fs | %s', ...
            char(method),k,toc(timer),ME.message));
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
        'CommitN:ResultAudit','Run result invalid for %s restart %d.',method,k);
    AMI_functions.assert_constellation_feasible(r.xBest,cfg,1e-10);
end

function path=find_m8(provided)
    if strlength(string(provided))>0
        path=char(provided);
        assert(isfile(path),'CommitN:M8Path','Provided M8ResultPath not found.');
        return;
    end
    root=fullfile(fso_result_utils.revision2_results_root(),'commit_m_final_fig7');
    hits=dir(fullfile(root,'**','CommitM_results.mat'));
    assert(~isempty(hits),'CommitN:MissingM8', ...
        'Please pass M8ResultPath pointing to your finished CommitM_results.mat.');
    [~,ii]=max([hits.datenum]);
    path=fullfile(hits(ii).folder,hits(ii).name);
end

function delete_tmp(file)
    if isfile(file)
        try,delete(file);catch,end
    end
end
function tf=is_text(x),tf=ischar(x)||(isstring(x)&&isscalar(x));end
function tf=is_bool(x),tf=islogical(x)&&isscalar(x);end
function tf=positive_int(x),tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x);end
