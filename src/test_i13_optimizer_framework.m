function report = test_i13_optimizer_framework()
%TEST_I13_OPTIMIZER_FRAMEWORK  Structural tests for Reviewer-2 framework.
%
% Does not require Global Optimization Toolbox; the actual Pattern Search
% execution is exercised by sim_i13_optimizer_pilot on the user's MATLAB.

    fprintf('\n============================================================\n');
    fprintf('I.13C - alternative optimizer framework regression\n');
    fprintf('============================================================\n');

    tests={@test_shared_starts,@test_start_projection,@test_fminsearch_smoke};
    names={ ...
        'Shared raw starts reproduce SA Threefry construction', ...
        'Shared starts project into the frozen feasible set', ...
        'Nelder-Mead projected alternative optimizer smoke test'};

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
        error('test_i13_optimizer_framework:Failure','I.13C framework regression failed.');
    end
end

function test_shared_starts()
    cfg=build_fso_config(8,1,20,0.2,'minGap',0.01,'pinZero',true, ...
        'saMaxIter',10,'saNStarts',6,'logEvery',0);
    base=linspace(0,2,8).';
    seed=uint32(1234567);
    X=i13_shared_start_points(cfg,base,seed,6);

    assert(isequal(X(:,1),base));
    s=RandStream('Threefry','Seed',double(seed));
    for k=2:6
        s.Substream=k;
        ref=base+0.2*randn(s,8,1);
        assert(max(abs(X(:,k)-ref))<1e-14);
    end
end

function test_start_projection()
    cfg=build_fso_config(32,1,30,0.2,'minGap',0.01,'pinZero',true, ...
        'saMaxIter',10,'saNStarts',8,'logEvery',0);
    base=linspace(0,2,32).';
    X=i13_shared_start_points(cfg,base,uint32(7654321),8);
    for k=1:size(X,2)
        xp=AMI_functions.project_constellation_1D(X(:,k),cfg);
        AMI_functions.assert_constellation_feasible(xp,cfg,1e-10);
    end
end

function test_fminsearch_smoke()
    cfg=build_fso_config(4,1,10,0, ...
        'minGap',0.01,'pinZero',true,'saMaxIter',10,'saNStarts',1,'logEvery',0);
    base=linspace(0,2,4).';
    [out,runs]=fminsearch_multistart_validated(cfg,base,uint32(24680),1,30);
    assert(isscalar(out.bestMIValidated)&&isfinite(out.bestMIValidated));
    assert(numel(runs)==1);
    AMI_functions.assert_constellation_feasible(out.bestXValidated,cfg,1e-8);
    assert(out.totalFunctionEvaluations<=35); % allow small MATLAB bookkeeping variation
end
