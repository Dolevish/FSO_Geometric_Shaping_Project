function report = test_i13d_final_optimizer_benchmark()
%TEST_I13D_FINAL_OPTIMIZER_BENCHMARK
% Structural regression for the frozen I.13D protocol.
%
% This test does not launch the expensive final benchmark. It verifies the
% pre-declared grid, restart/evaluation budgets, unique task identities and
% longest-first priority policy from a synthetic full I.11.1 case table.

    fprintf('\n============================================================\n');
    fprintf('I.13D - final optimizer benchmark structural regression\n');
    fprintf('============================================================\n');

    tests={@test_frozen_config,@test_task_plan,@test_priority_policy,@test_tiny_final_repair};
    names={ ...
        'Frozen case/restart/evaluation counts', ...
        '32-case / 224-task plan with unique identities', ...
        'Longest-first plan prioritizes turbulent refined high-M work', ...
        'Tiny Pattern Search d_min residual is canonically repaired'};

    passed=false(numel(tests),1); elapsed=nan(numel(tests),1); messages=strings(numel(tests),1);
    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k}); t=tic;
        try
            tests{k}(); passed(k)=true; messages(k)="PASS";
            elapsed(k)=toc(t); fprintf('PASS (%.2fs)\n',elapsed(k));
        catch ME
            elapsed(k)=toc(t); messages(k)=string(ME.message);
            fprintf('FAIL (%.2fs)\n  %s\n',elapsed(k),ME.message);
        end
    end

    report=struct('names',{names},'passed',passed,'elapsedSeconds',elapsed, ...
        'messages',messages,'nPassed',sum(passed),'nTotal',numel(tests));
    fprintf('\nPassed %d/%d tests in %.2fs.\n',report.nPassed,report.nTotal,sum(elapsed,'omitnan'));
    if ~all(passed)
        error('test_i13d_final_optimizer_benchmark:Failure','I.13D structural regression failed.');
    end
end


function test_frozen_config()
    F=i13d_benchmark_config();
    assert(isequal(F.MVec,[4 8 16 32]));
    assert(isequal(F.SNRVec,[20 30]));
    assert(isequal(F.TurbulenceVec,[0 0.3]));
    assert(isequal(F.Replicates,[1 2]));
    assert(isequal(F.RestartsByM,[6 6 8 8]));
    assert(isequal(F.MaxIterByM,[4000 4000 4000 8000]));
    assert(isequal(F.MaxEvalsPerStart,[4001 4001 4001 8001]));
    assert(F.nPhysicalPoints==16);
    assert(F.nOptimizationCases==32);
    assert(F.nRestartTasks==224);
    assert(F.totalMaxFunctionEvaluations==1152224);
    assert(abs(F.tieToleranceBits-1e-3)<eps);
    assert(strcmp(F.scheduler,'parfeval-dynamic-v1'));
    assert(~F.PatternSearch.UseParallel);
    assert(~F.PatternSearch.UseCompletePoll);
end


function test_task_plan()
    F=i13d_benchmark_config();
    T=synthetic_case_table();
    [C,Q]=i13d_build_task_plan(T,F);

    assert(height(C)==32);
    assert(height(Q)==224);
    assert(numel(unique(Q.IdentityKey))==224);

    for M=F.MVec
        iM=find(F.MVec==M,1);
        c=C.M==M;
        assert(all(C.NStarts(c)==F.RestartsByM(iM)));
        assert(all(C.MaxEvalsPerStart(c)==F.MaxEvalsPerStart(iM)));
    end
end


function test_priority_policy()
    F=i13d_benchmark_config();
    [~,Q]=i13d_build_task_plan(synthetic_case_table(),F);

    key=char(Q.IdentityKey(1));
    assert(contains(key,'M32_SNR30_sig300'));
    assert(all(diff(Q.Priority)<=0));
end


function test_tiny_final_repair()
    cfg=build_fso_config(4,1,20,0.3, ...
        'minGap',0.01,'pinZero',true,'saMaxIter',1,'saNStarts',1,'logEvery',0);

    % Exact feasible boundary constellation: mean=1, x1=0, first two
    % adjacent gaps exactly d_min.
    x=[0;0.01;0.02;3.97];

    % Mimic the observed Pattern Search numerical residual while preserving
    % mean power: first gap is short by 9.7e-7.
    xBad=x;
    xBad(2)=xBad(2)-9.7e-7;
    xBad(4)=xBad(4)+9.7e-7;

    [xFix,meta]=i13d_finalize_patternsearch_solution(xBad,cfg, ...
        'MaxRepairInfNorm',1e-4);

    AMI_functions.assert_constellation_feasible(xFix,cfg,1e-10);
    assert(meta.wasRepaired);
    assert(meta.before.minGapDeficit>9e-7 && meta.before.minGapDeficit<1.1e-6);
    assert(meta.repairInfNorm<1e-4);
end


function T=synthetic_case_table()
    P=i11p1_production_config();
    n=P.nOptimizationCases;

    M=zeros(n,1); SigmaX2=zeros(n,1); SNRdB=zeros(n,1); Replicate=zeros(n,1);
    Seed=zeros(n,1); NStarts=zeros(n,1); SelectedDt=zeros(n,1);
    AMIValidated=zeros(n,1); GainBits=zeros(n,1);
    q=0;

    for m=P.MVec
        [starts,~]=i11_restart_policy(m);
        for sig=P.TurbulenceVec
            for snr=P.SNRVec
                for rep=1:P.nReplicates
                    q=q+1;
                    M(q)=m; SigmaX2(q)=sig; SNRdB(q)=snr; Replicate(q)=rep;
                    Seed(q)=100000+q; NStarts(q)=starts;
                    hard = snr>=30 && ((m>=32 && sig>=0.1) || (m>=16 && sig>=0.2));
                    if hard,SelectedDt(q)=0.005;else,SelectedDt(q)=0.01;end
                    AMIValidated(q)=0.5+0.01*m+0.001*snr-0.1*sig+0.001*rep;
                    GainBits(q)=0.1;
                end
            end
        end
    end

    T=table(M,SigmaX2,SNRdB,Replicate,Seed,NStarts,SelectedDt,AMIValidated,GainBits);
end
