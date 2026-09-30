function report = test_i8_high_snr_accuracy()
%TEST_I8_HIGH_SNR_ACCURACY  Regression checks for Commit I.8 diagnostic.

    fprintf('\n============================================================\n');
    fprintf('Commit I.8 - high-SNR accuracy diagnostic regression suite\n');
    fprintf('============================================================\n');

    tests={@test_smoke_factorial,@test_production_defaults_untouched};
    names={'Factorial dt/y refinement smoke','I.8 does not alter production defaults'};
    passed=false(numel(tests),1); elapsed=nan(numel(tests),1); messages=strings(numel(tests),1);

    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k});
        t=tic;
        try
            tests{k}();
            elapsed(k)=toc(t); passed(k)=true; messages(k)="PASS";
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
        error('test_i8_high_snr_accuracy:Failure','Commit I.8 regression suite failed.');
    end
end


function test_smoke_factorial()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    d=diagnose_high_snr_logh_accuracy( ...
        'Cases',[4 20 0.1], ...
        'GeometryMode','uniform', ...
        'DtVec',[0.01 0.005], ...
        'YRefineVec',[1 2], ...
        'UseParallelCases',false, ...
        'WriteCSV',false, ...
        'ResultsRoot',root, ...
        'RunLabel','I8_regression');

    T=d.results;
    assert(height(T)==8,'Expected 2 modes x 2 dt x 2 y-refinement rows.');
    assert(all(T.AbsError<1e-3),'Representative smoke case should satisfy accuracy target.');
    assert(all(abs(T.ValidatorAMI-T.ValidatorAMI(1))<1e-14),'All settings must share one validator value.');

    for mode=["fixed-bound","candidate-adaptive"]
        b=T(T.Mode==mode & T.YRefine==1 & abs(T.Dt-0.01)<1e-15,:);
        r=T(T.Mode==mode & T.YRefine==2 & abs(T.Dt-0.01)<1e-15,:);
        assert(height(b)==1 && height(r)==1);
        assert(abs(b.YSupportMax-r.YSupportMax)<1e-12,'y refinement must preserve support endpoint.');
        assert(r.YGridN==2*(b.YGridN-1)+1,'factor-2 refinement must only subdivide intervals.');
    end

    assert(d.summary.coarsestPassingAdaptiveDtAtBaseY>=0.01-1e-15);
end


function test_production_defaults_untouched()
    % I.8 is diagnostic-only. I.7 production settings must remain explicit.
    cfg=build_fso_config(4,1,20,0.1,'minGap',0.01,'saMaxIter',1,'saNStarts',1);
    assert(strcmp(cfg.yGridMode,'fixed-bound'), ...
        'Canonical build_fso_config default must remain fixed-bound for reproducibility.');
    assert(abs(cfg.loghDt-0.01)<1e-14,'Canonical log-h dt must remain 0.01 in I.8.');
end


function safe_rmdir(p)
    if exist(p,'dir')
        try, rmdir(p,'s'); catch, end
    end
end
