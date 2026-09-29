function report = test_i7_production_integration()
%TEST_I7_PRODUCTION_INTEGRATION  Regression checks for Commit I.7.
%
% Verifies that the production AMI driver explicitly selects the validated
% candidate-adaptive y-grid, that metadata/case snapshots retain this choice,
% and that the no-SA production-grid accuracy diagnostic executes and meets
% the configured tolerance on a representative smoke point.

    fprintf('\n============================================================\n');
    fprintf('Commit I.7 - production integration regression suite\n');
    fprintf('============================================================\n');

    tests = {@test_production_driver_default,@test_accuracy_sweep_smoke, ...
        @test_fixed_reproducibility_path,@test_gh_adaptive_rejected};
    names = {'Production driver selects adaptive grid', ...
        'Full-grid diagnostic smoke passes tolerance', ...
        'Fixed-bound reproducibility path remains available', ...
        'GH plus adaptive grid is rejected explicitly'};

    passed = false(numel(tests),1);
    elapsed = nan(numel(tests),1);
    messages = strings(numel(tests),1);

    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k});
        t=tic;
        try
            tests{k}();
            elapsed(k)=toc(t); passed(k)=true; messages(k)="PASS";
            fprintf('PASS (%.2fs)\n',elapsed(k));
        catch ME
            elapsed(k)=toc(t); messages(k)=string(ME.message);
            fprintf('FAIL (%.2fs)\n',elapsed(k));
            fprintf('  %s\n',ME.message);
        end
    end

    report = struct('names',{names},'passed',passed,'elapsedSeconds',elapsed, ...
        'messages',messages,'nPassed',sum(passed),'nTotal',numel(tests));
    fprintf('\nPassed %d/%d tests in %.2fs.\n',report.nPassed,report.nTotal,sum(elapsed,'omitnan'));
    if ~all(passed)
        error('test_i7_production_integration:Failure','Commit I.7 regression suite failed.');
    end
end


function test_production_driver_default()
    root=tempname;
    mkdir(root);
    cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>

    b=sim_revision_AMI_vs_SNR( ...
        'MVec',4,'SNRVec',20,'TurbulenceVec',0.1, ...
        'saMaxIter',1,'saNStarts',1,'nReplicates',1, ...
        'UseParallelCases',false,'MakePlots',false,'SaveFigures',false, ...
        'ResultsRoot',root,'RunLabel','I7_regression');

    assert(strcmp(b.method.fastFadingMethod,'logh'));
    assert(strcmp(b.method.yGridMode,'candidate-adaptive'));
    assert(abs(b.method.yGridScale-1.1)<1e-14);
    assert(contains(string(b.run.methodTag),'Yadapt1p10'));
    assert(abs(b.summary.pamFastValidationGap(1))<1e-3);

    p=b.caseFiles{1,1,1,1};
    s=load(p,'c');
    assert(strcmp(s.c.config.yGridMode,'candidate-adaptive'));
    assert(abs(s.c.config.yGridScale-1.1)<1e-14);
end


function test_accuracy_sweep_smoke()
    root=tempname;
    mkdir(root);
    cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>

    d=diagnose_production_y_grid_accuracy( ...
        'MVec',4,'SNRVec',20,'TurbulenceVec',0.1, ...
        'GeometryModes',{'uniform','upper-eighth'}, ...
        'UseParallelCases',false,'WriteCSV',false, ...
        'ResultsRoot',root,'RunLabel','I7_diag_regression');

    assert(d.summary.adaptiveAllPass);
    assert(d.summary.nAdaptiveFailures==0);
    assert(d.summary.maxAbsAdaptiveError<=d.method.absErrorTolerance);
    assert(strcmp(d.method.adaptiveYGridMode,'candidate-adaptive'));
    assert(abs(d.method.adaptiveScale-1.1)<1e-14);
end


function test_fixed_reproducibility_path()
    cfg=build_fso_config(8,1,20,0.1, ...
        'FastFadingMethod','logh','YGridMode','fixed-bound','YGridScale',1.1, ...
        'xMaxBound',fso_result_utils.feasible_xmax_bound(8,1,0.01,true), ...
        'minGap',0.01,'pinZero',true,'saMaxIter',1,'saNStarts',1);
    assert(strcmp(cfg.yGridMode,'fixed-bound'));
    x=define_constellation(8,1,true,"mean");
    a=cfg.AMI_Evaluator(x);
    b=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,cfg.loghDt,cfg.loghSpanSigma, ...
        cfg.y_grid,cfg.loghYBlockSize);
    assert(abs(a-b)<1e-12);
end


function test_gh_adaptive_rejected()
    threw=false;
    try
        build_fso_config(4,1,20,0.1, ...
            'FastFadingMethod','gh','ghN_h',40, ...
            'YGridMode','candidate-adaptive','minGap',0.01);
    catch ME
        threw=strcmp(ME.identifier,'build_fso_config:AdaptiveYGridRequiresLogh');
    end
    assert(threw,'Expected explicit GH/adaptive incompatibility rejection.');
end


function safe_rmdir(p)
    if exist(p,'dir')
        try
            rmdir(p,'s');
        catch
        end
    end
end
