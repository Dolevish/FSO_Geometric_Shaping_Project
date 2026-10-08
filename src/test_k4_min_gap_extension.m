function report = test_k4_min_gap_extension()
%TEST_K4_MIN_GAP_EXTENSION  Structural/synthetic preflight, no SA run.
%
% Validates the six feasible points, paired seeds, unchanged base protocol,
% dry-run task dispatch, and synthetic 160+96=256 restart table merge.
% Tests combined MATLAB figure exports without running the AMI evaluators.

    fprintf('\n============================================================\n');
    fprintf('COMMIT K4 - extension preflight\n');
    fprintf('============================================================\n');

    tests={@test_config,@test_plan,@test_pairing,@test_original_unchanged, ...
        @test_dry_run,@test_merge_and_figure};
    names={'Six feasible points and compute budget', ...
        '12 cases / 96 unique restarts; no infeasible M32 points', ...
        'Paired case seeds match the completed K protocol', ...
        'Original K-v1 grid, counts, and task identities unchanged', ...
        'K4 engine dry-run uses the six-point plan', ...
        'Synthetic table merge and publication-figure export'};

    passed=false(numel(tests),1);
    messages=strings(numel(tests),1);
    elapsed=zeros(numel(tests),1);

    for k=1:numel(tests)
        fprintf('[%d/%d] %s ... ',k,numel(tests),names{k});
        t=tic;
        try
            tests{k}();
            passed(k)=true;
            messages(k)="PASS";
            fprintf('PASS (%.2f s)\n',toc(t));
        catch ME
            messages(k)=string(ME.message);
            fprintf('FAIL (%s)\n',ME.message);
        end
        elapsed(k)=toc(t);
    end

    report=struct('passed',passed,'messages',messages,'elapsedSeconds',elapsed, ...
        'nPassed',nnz(passed),'nTotal',numel(tests));
    fprintf('\nCommit K4 preflight: %d/%d PASS\n',report.nPassed,report.nTotal);
    if ~all(passed)
        error('test_k4_min_gap_extension:Failed', ...
            '%d of %d K4 preflight tests failed.',nnz(~passed),numel(tests));
    end
end


function test_config()
    K=k4_min_gap_extension_config();
    assert(isequal(K.MVec,[16 32]));
    assert(isequal(K.MinGapVec,[0.03 0.04 0.07 0.10]));
    assert(isequal(K.MinGapByM{1},[0.03 0.04 0.07 0.10]));
    assert(isequal(K.MinGapByM{2},[0.03 0.04]));
    assert(isequal(K.Replicates,[1 2]));
    assert(K.NStarts==8 && K.MaxIter==8000);
    assert(K.nPhysicalSensitivityPoints==6);
    assert(K.nReplicateCases==12);
    assert(K.nRestartTasks==96);
    assert(K.totalSAIterations==768000);
    assert(K.nominalFastObjectiveEvaluations==768096);
    assert(all(abs(K.SelectedDtByM-0.01)<1e-12));
    assert(0.07>2/(32-1) && 0.10>2/(32-1));
end


function test_plan()
    K=k4_min_gap_extension_config();
    [C,T]=k_build_task_plan(K);
    assert(height(C)==12 && height(T)==96);
    assert(numel(unique(T.IdentityKey))==96);
    assert(all(C.MinGap(C.M==32)<=2/(32-1)+1e-12));
    assert(~any(C.M==32 & C.MinGap>=0.07));
    assert(nnz(C.M==16)==8 && nnz(C.M==32)==4);
    first=C(T.TaskCaseIndex(1),:);
    assert(first.M==32 && abs(first.MinGap-0.03)<1e-12);

    for c=1:height(C)
        assert(nnz(T.TaskCaseIndex==c)==8);
        assert(all(sort(T.Restart(T.TaskCaseIndex==c)).'==1:8));
        x=define_constellation(C.M(c),K.P_avg,true,"mean");
        assert(min(diff(x))>=C.MinGap(c)-1e-12);
    end
end


function test_pairing()
    K=k_min_gap_config();
    E=k4_min_gap_extension_config();
    [C,~]=k_build_task_plan(K);
    [D,~]=k_build_task_plan(E);
    for m=K.MVec
        for rep=K.Replicates
            a=unique(C.Seed(C.M==m & C.Replicate==rep));
            b=unique(D.Seed(D.M==m & D.Replicate==rep));
            assert(numel(a)==1 && numel(b)==1 && a==b);
        end
    end
end


function test_original_unchanged()
    K=k_min_gap_config();
    assert(strcmp(K.version,'K-v1'));
    assert(~isfield(K,'MinGapByM'));
    assert(isequal(K.MinGapVec,[0 0.005 0.01 0.02 0.05]));
    assert(K.nPhysicalSensitivityPoints==10 && K.nRestartTasks==160);
    [C,T]=k_build_task_plan(K);
    assert(height(C)==20 && height(T)==160);
    assert(numel(unique(T.IdentityKey))==160);
end


function test_dry_run()
    txt=evalc("R4=sim_k4_min_gap_extension('DryRun',true);");
    assert(R4.dryRun);
    assert(strcmp(R4.version,'K4-v1'));
    assert(height(R4.casePlan)==12 && height(R4.taskPlan)==96);
    assert(contains(txt,'restart tasks=96'));
end


function test_merge_and_figure()
    K=k_min_gap_config();
    E=k4_min_gap_extension_config();
    A=synthetic_result(K);
    B=synthetic_result(E);

    resultsRoot=fso_result_utils.revision2_results_root();
    if exist(resultsRoot,'dir')~=7,mkdir(resultsRoot);end
    tempDir=tempname(resultsRoot);
    mkdir(tempDir);
    cleanup=onCleanup(@() cleanup_dir(tempDir)); %#ok<NASGU>

    C=merge_k4_min_gap_results(A,B,'OutputDirectory',tempDir);
    assert(height(C.physicalSummary)==16);
    assert(height(C.replicateSummary)==32);
    assert(height(C.restartSummary)==256);
    assert(nnz(C.physicalSummary.M==16)==9);
    assert(nnz(C.physicalSummary.M==32)==7);
    assert(all(isfinite(C.physicalSummary.GainRetentionVsGap0)));
    assert(exist(fullfile(tempDir,'K4_combined_results.mat'),'file')==2);
    assert(exist(fullfile(tempDir,'K4_combined_physical_summary.csv'),'file')==2);
    assert(exist(C.figurePaths.pdf,'file')==2);
    assert(exist(C.figurePaths.png,'file')==2);
    assert(exist(C.figurePaths.fig,'file')==2);
    close all;
end


function S=synthetic_result(K)
    [cases,tasks]=k_build_task_plan(K);
    n=height(tasks);
    M=zeros(n,1);MinGap=zeros(n,1);Replicate=zeros(n,1);
    Restart=zeros(n,1);ValidatedAMI=zeros(n,1);FastValidatorGap=zeros(n,1);
    for i=1:n
        C=cases(tasks.TaskCaseIndex(i),:);
        M(i)=C.M;MinGap(i)=C.MinGap;Replicate(i)=C.Replicate;
        Restart(i)=tasks.Restart(i);
        ValidatedAMI(i)=synthetic_ami(C.M,C.MinGap,C.Replicate) - ...
            0.00001*(8-Restart(i));
    end
    restarts=table(M,MinGap,Replicate,Restart,ValidatedAMI,FastValidatorGap);
    restarts=sortrows(restarts,{'M','MinGap','Replicate','Restart'});

    r=height(cases);
    M=zeros(r,1);MinGap=zeros(r,1);Replicate=zeros(r,1);
    WinnerRestart=8*ones(r,1);WinnerValidatedAMI=zeros(r,1);
    WinnerConstellation=cell(r,1);SelectionChanged=false(r,1);
    for k=1:r
        C=cases(k,:);
        M(k)=C.M;MinGap(k)=C.MinGap;Replicate(k)=C.Replicate;
        WinnerValidatedAMI(k)=synthetic_ami(C.M,C.MinGap,C.Replicate);
        WinnerConstellation{k}=define_constellation(C.M,K.P_avg,true,"mean");
    end
    replicateSummary=table(M,MinGap,Replicate,WinnerRestart, ...
        WinnerValidatedAMI,SelectionChanged,WinnerConstellation);
    replicateSummary=sortrows(replicateSummary,{'M','MinGap','Replicate'});

    keys=unique(replicateSummary(:,{'M','MinGap'}),'rows','sorted');
    n=height(keys);
    Rep1ValidatedAMI=zeros(n,1);Rep2ValidatedAMI=zeros(n,1);
    MeanValidatedAMI=zeros(n,1);StdAcrossReplicates=zeros(n,1);
    ReplicateRange=zeros(n,1);MeanGainOverPAM=zeros(n,1);
    GainRetentionVsGap0=nan(n,1);
    for k=1:n
        m=keys.M(k);d=keys.MinGap(k);
        Rep1ValidatedAMI(k)=synthetic_ami(m,d,1);
        Rep2ValidatedAMI(k)=synthetic_ami(m,d,2);
        MeanValidatedAMI(k)=mean([Rep1ValidatedAMI(k) Rep2ValidatedAMI(k)]);
        StdAcrossReplicates(k)=std([Rep1ValidatedAMI(k) Rep2ValidatedAMI(k)],0);
        ReplicateRange(k)=abs(Rep1ValidatedAMI(k)-Rep2ValidatedAMI(k));
        MeanGainOverPAM(k)=MeanValidatedAMI(k)-0.8;
        if strcmp(K.version,'K-v1')
            GainRetentionVsGap0(k)=MeanGainOverPAM(k)/ ...
                (mean([synthetic_ami(m,0,1) synthetic_ami(m,0,2)])-0.8);
        end
    end
    physicalSummary=[keys table(Rep1ValidatedAMI,Rep2ValidatedAMI, ...
        MeanValidatedAMI,StdAcrossReplicates,ReplicateRange, ...
        MeanGainOverPAM,GainRetentionVsGap0)];

    B=struct('M',cell(2,1),'x',[],'amiFast',[],'amiValidated',[]);
    for i=1:2
        m=K.MVec(i);
        B(i).M=m;
        B(i).x=define_constellation(m,K.P_avg,true,"mean");
        B(i).amiFast=0.8;
        B(i).amiValidated=0.8;
    end

    S=struct();
    S.version=K.version; S.config=K;S.baselines=B;
    S.casePlan=cases;S.restartSummary=restarts;
    S.replicateSummary=replicateSummary;S.physicalSummary=physicalSummary;
    S.metrics=struct('scheduler',struct('workerUtilization',1));
    S.outputDirectory=fso_result_utils.revision2_results_root();
end


function v=synthetic_ami(m,d,rep)
    v=1.30+0.003*m-0.9*d+0.0005*rep;
end


function cleanup_dir(path)
    if exist(path,'dir')==7
        try,rmdir(path,'s');catch,end
    end
end
