function report = test_i13_zhou_baseline()
%TEST_I13_ZHOU_BASELINE  Regression tests for the Reviewer-1 GS baseline.

    fprintf('\n============================================================\n');
    fprintf('I.13B - Zhou-2025 literature GS baseline regression\n');
    fprintf('============================================================\n');

    tests={@test_integer_levels,@test_power_and_feasibility,@test_validator_smoke};
    names={ ...
        'Published-equation n=2 integer levels for M=4,8,16,32', ...
        'Mean power, zero anchor, ordering and frozen d_min feasibility', ...
        'Independent FSO validator smoke test'};

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
        error('test_i13_zhou_baseline:Failure','I.13B baseline regression failed.');
    end
end

function test_integer_levels()
    expected=cell(4,1);
    expected{1}=[0 2 6 15];
    expected{2}=[0 1 3 5 8 11 17 31];
    expected{3}=[0 1 2 4 5 7 8 10 12 15 17 21 25 31 40 63];
    expected{4}=[0 1 2 3 4 5 6 7 8 10 11 12 14 15 17 18 20 22 24 26 29 31 34 37 41 45 50 56 63 73 87 127];
    Mvec=[4 8 16 32];

    for k=1:numel(Mvec)
        Z=zhou2025_optical_gs_constellation(Mvec(k),1,'ExtraBits',2);
        assert(isequal(Z.integerLevels(:).',expected{k}), ...
            'Unexpected n=2 integer levels for M=%d.',Mvec(k));
        assert(Z.integerLevels(end)==2^(log2(Mvec(k))+2)-1);
    end
end

function test_power_and_feasibility()
    Mvec=[4 8 16 32];
    for M=Mvec
        Z=zhou2025_optical_gs_constellation(M,1,'ExtraBits',2);
        assert(abs(mean(Z.x)-1)<1e-12);
        assert(abs(Z.x(1))<1e-12);
        assert(all(diff(Z.x)>0));
        assert(Z.minGap>=0.01);

        cfg=build_fso_config(M,1,20,0.2, ...
            'minGap',0.01,'pinZero',true,'saMaxIter',1,'saNStarts',1,'logEvery',0);
        AMI_functions.assert_constellation_feasible(Z.x,cfg,1e-10);
    end
end

function test_validator_smoke()
    Z=zhou2025_optical_gs_constellation(4,1,'ExtraBits',2);
    cfg=build_fso_config(4,1,15,0, ...
        'minGap',0.01,'pinZero',true,'saMaxIter',1,'saNStarts',1,'logEvery',0);
    I=cfg.AMI_Validator(Z.x);
    assert(isscalar(I)&&isfinite(I)&&I>=0&&I<=log2(4)+1e-6);
end
