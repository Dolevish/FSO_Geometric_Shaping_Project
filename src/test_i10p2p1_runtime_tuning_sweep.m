function report = test_i10p2p1_runtime_tuning_sweep()
%TEST_I10P2P1_RUNTIME_TUNING_SWEEP  Regression checks for tuning sweep/freeze.

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.2.1/I.10.2.2 - runtime tuning regression suite\n');
    fprintf('============================================================\n');

    tests={@test_selector_matches_frozen_policy,@test_small_sweep_assembles};
    names={'Production runtime selector matches frozen measured policy', ...
        'Small tuning sweep builds robust summaries without changing AMI'};

    passed=false(numel(tests),1); elapsed=nan(numel(tests),1); messages=strings(numel(tests),1);
    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k}); t=tic;
        try
            tests{k}(); elapsed(k)=toc(t); passed(k)=true; messages(k)="PASS";
            fprintf('PASS (%.2fs)\n',elapsed(k));
        catch ME
            elapsed(k)=toc(t); messages(k)=string(ME.message);
            fprintf('FAIL (%.2fs)\n  %s\n',elapsed(k),ME.message);
        end
    end

    report=struct('names',{names},'passed',passed,'elapsedSeconds',elapsed, ...
        'messages',messages,'nPassed',sum(passed),'nTotal',numel(tests));
    fprintf('\nPassed %d/%d tests in %.2fs.\n',report.nPassed,report.nTotal,sum(elapsed,'omitnan'));
    if ~all(passed)
        error('test_i10p2p1_runtime_tuning_sweep:Failure', ...
            'Commit I.10.2.1/I.10.2.2 runtime tuning regression suite failed.');
    end
end


function test_selector_matches_frozen_policy()
    expected=[4 2048 1;8 512 4;16 256 4;32 256 4];
    for k=1:size(expected,1)
        M=expected(k,1); [b,s,m]=select_logh_runtime_tuning(M);
        assert(b==expected(k,2) && s==expected(k,3), ...
            'Unexpected frozen runtime setting for M=%d.',M);
        assert(strcmp(m.ruleVersion,'I10.2.2-v1'), ...
            'Unexpected runtime tuning rule version for M=%d.',M);
        assert(abs(m.noRegressionThreshold-0.98)<1e-12, ...
            'Unexpected no-regression threshold metadata for M=%d.',M);
    end
end


function test_small_sweep_assembles()
    root=fullfile(tempdir,'fso_i10p2p1_regression');
    d=diagnose_i10p2p1_runtime_tuning_sweep( ...
        'Cases',[4 20 0.1 0.01;32 30 0.3 0.005], ...
        'GeometryModes',"uniform", ...
        'BlockSizes',[256 512], ...
        'BatchSizes',[1 2], ...
        'nWarmup',0, ...
        'nRepeats',1, ...
        'WriteCSV',false, ...
        'ResultsRoot',root, ...
        'RunLabel','regression');

    assert(d.summary.nPhysicalCases==2,'Expected two physical cases.');
    assert(d.summary.nScenarios==2,'Expected one geometry per physical case.');
    assert(height(d.results)==10, ...
        'Expected auto plus four tuned rows for each of two scenarios.');
    assert(d.summary.maxAbsDifference<=5e-12, ...
        'Tuning sweep changed AMI beyond numerical tolerance.');
    assert(height(d.summary.recommendedByM)==2,'Expected recommendation rows for M=4 and M=32.');
    assert(height(d.summary.bestMedianByM)==2,'Expected best-median rows for M=4 and M=32.');
    assert(all(ismember(d.summary.recommendedByM.M,[4;32])), ...
        'Unexpected M values in recommendation table.');
end
