function report = test_i10_final_sa_pilot()
%TEST_I10_FINAL_SA_PILOT  Regression checks for Commit I.10.

    fprintf('\n============================================================\n');
    fprintf('Commit I.10 - final SA pilot regression suite\n');
    fprintf('============================================================\n');

    tests={@test_default_case_policy_map,@test_boundary_pair_smoke};
    names={'Default five-case policy map is exact', ...
        'Boundary paired smoke uses identical seeds and dt 0.01 vs 0.005'};

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
    if ~all(passed), error('test_i10_final_sa_pilot:Failure','Commit I.10 regression suite failed.'); end
end


function test_default_case_policy_map()
    C=[32 30 0.1;32 30 0.2;32 30 0.3;16 30 0.3;16 30 0.2];
    expected=[0.005;0.005;0.005;0.005;0.01];
    got=nan(size(expected)); refined=false(size(expected));
    for k=1:size(C,1)
        [got(k),m]=select_logh_dt(C(k,1),C(k,2),C(k,3)); refined(k)=m.isRefined;
    end
    assert(max(abs(got-expected))<1e-14,'Unexpected dt mapping for the five I.10 cases.');
    assert(isequal(refined,[true;true;true;true;false]),'Unexpected refined mask.');
end


function test_boundary_pair_smoke()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>

    % Low-complexity structural smoke: same A/B mechanism as the real boundary
    % comparison, but on M=4 to keep the regression test fast.
    b=sim_i10_final_sa_pilot( ...
        'Cases',[4 20 0.1], ...
        'BoundaryCase',[4 20 0.1], ...
        'RunBoundaryReference',true, ...
        'saMaxIter',1,'saNStarts',2,'nReplicates',1, ...
        'UseParallelChildRuns',false,'WriteCSV',false, ...
        'ResultsRoot',root,'RunLabel','I10_regression');

    T=b.results; P=b.boundaryPairs;
    assert(height(T)==2,'Expected one policy row and one boundary-reference row.');
    assert(height(P)==1,'Expected exactly one paired boundary comparison.');
    a=T(T.Mode=="policy",:); r=T(T.Mode=="boundary-reference",:);
    assert(height(a)==1 && height(r)==1);
    assert(abs(a.SelectedDt-0.01)<1e-14,'Smoke policy case should use base dt=0.01.');
    assert(abs(r.SelectedDt-0.005)<1e-14,'Reference must use fixed dt=0.005.');
    assert(a.Seed==r.Seed && P.SameSeed,'Policy/reference runs must use identical seeds.');
    assert(strcmp(b.policy.ruleVersion,'I9-v1'));
end


function safe_rmdir(p)
    if exist(p,'dir'), try, rmdir(p,'s'); catch, end; end
end
