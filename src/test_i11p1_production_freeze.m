function report = test_i11p1_production_freeze()
%TEST_I11P1_PRODUCTION_FREEZE  Regression checks for the revised production freeze.

    fprintf('\n============================================================\n');
    fprintf('Commit I.11.1 - variable-iteration production freeze suite\n');
    fprintf('============================================================\n');

    tests={@test_iteration_policy,@test_frozen_profile,@test_stage_plan_and_alias};
    names={'Per-M iteration policy is exactly 4000/4000/4000/8000', ...
        'Frozen production totals and restart map are exact', ...
        'Canonical stage plan runs M32 first and legacy alias targets I.11.1'};

    passed=false(numel(tests),1); elapsed=nan(numel(tests),1); messages=strings(numel(tests),1);
    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k}); t=tic;
        try
            tests{k}(); passed(k)=true; messages(k)="PASS";
            fprintf('PASS (%.2fs)\n',toc(t));
        catch ME
            elapsed(k)=toc(t); messages(k)=string(ME.message);
            fprintf('FAIL (%.2fs)\n  %s\n',elapsed(k),ME.message);
            continue;
        end
        elapsed(k)=toc(t);
    end

    report=struct('names',{names},'passed',passed,'elapsedSeconds',elapsed, ...
        'messages',messages,'nPassed',sum(passed),'nTotal',numel(tests));
    fprintf('\nPassed %d/%d tests in %.2fs.\n',report.nPassed,report.nTotal,sum(elapsed,'omitnan'));
    if ~all(passed)
        error('test_i11p1_production_freeze:Failure','I.11.1 production-freeze regression suite failed.');
    end
end

function test_iteration_policy()
    M=[4 8 16 32]; got=zeros(size(M));
    for k=1:numel(M)
        [got(k),meta]=i11p1_iteration_policy(M(k));
        assert(strcmp(meta.ruleVersion,'I11.1-v1'));
        assert(~meta.globalOptimumClaim);
    end
    assert(isequal(got,[4000 4000 4000 8000]));
end

function test_frozen_profile()
    F=i11p1_production_config();
    assert(strcmp(F.freezeVersion,'I11.1-v1'));
    assert(isequal(F.MVec,[4 8 16 32]));
    assert(isequal(F.RestartsByM,[6 6 8 8]));
    assert(isequal(F.MaxIterByM,[4000 4000 4000 8000]));
    assert(F.nReplicates==2);
    assert(F.nOptimizationCases==192);
    assert(F.totalRestartTasks==1344);
    assert(F.totalSAIterations==6912000);
    assert(strcmp(F.dtRuleVersion,'I10.3-v1'));
    assert(strcmp(F.runtimeTuningVersion,'I10.2.2-v1'));
    assert(~F.globalOptimumClaim);
end

function test_stage_plan_and_alias()
    F=i11p1_production_config(); P=i11p1_stage_plan(F);
    assert(numel(P)==2);
    assert(strcmp(P(1).name,'m32') && isequal(P(1).MVec,32) && P(1).saMaxIter==8000);
    assert(strcmp(P(2).name,'standard') && isequal(P(2).MVec,[4 8 16]) && P(2).saMaxIter==4000);
    assert(P(1).nOptimizationCases==48 && P(1).nRestartTasks==384);
    assert(P(2).nOptimizationCases==144 && P(2).nRestartTasks==960);
    assert(P(1).totalSAIterations==3072000);
    assert(P(2).totalSAIterations==3840000);

    b=sim_i11_production_freeze('DryRun',true);
    assert(isfield(b,'dryRun') && b.dryRun);
    assert(isequal(b.productionFreeze.MaxIterByM,[4000 4000 4000 8000]));
    assert(strcmp(b.stagePlan(1).name,'m32'));
end
