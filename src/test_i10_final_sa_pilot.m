function report = test_i10_final_sa_pilot()
%TEST_I10_FINAL_SA_PILOT  Regression checks for Commit I.10 / I.10.3 policy.

    fprintf('\n============================================================\n');
    fprintf('Commit I.10 / I.10.3 - final SA pilot regression suite\n');
    fprintf('============================================================\n');

    tests={@test_default_case_policy_map,@test_boundary_pair_smoke};
    names={'Default five-case revised policy map is exact', ...
        'Historical paired-smoke mechanism still uses identical seeds'};

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
    if ~all(passed), error('test_i10_final_sa_pilot:Failure','Commit I.10/I.10.3 regression suite failed.'); end
end


function test_default_case_policy_map()
    C=[32 30 0.1;32 30 0.2;32 30 0.3;16 30 0.3;16 30 0.2];
    expected=0.005*ones(5,1);
    got=nan(size(expected)); refined=false(size(expected)); versions=strings(size(expected));
    for k=1:size(C,1)
        [got(k),m]=select_logh_dt(C(k,1),C(k,2),C(k,3));
        refined(k)=m.isRefined; versions(k)=string(m.ruleVersion);
    end
    assert(max(abs(got-expected))<1e-14,'Unexpected dt mapping for the five revised policy cases.');
    assert(all(refined),'All five I.10/I.10.3 cases must now use refined dt=0.005.');
    assert(all(versions=="I10.3-v1"),'Unexpected revised dt-policy rule version.');
end


function test_boundary_pair_smoke()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>

    % Historical structural smoke for the paired A/B machinery.  We keep a
    % low-complexity M=4 base-policy case so policy dt remains 0.01 while the
    % explicit reference remains 0.005.  This validates the mechanism without
    % re-running the obsolete M16 boundary comparison after I.10.3 refined it.
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
    assert(strcmp(b.policy.selectedMeta{1}.ruleVersion,'I10.3-v1'), ...
        'Pilot must carry the revised selector metadata from select_logh_dt.');
end


function safe_rmdir(p)
    if exist(p,'dir'), try, rmdir(p,'s'); catch, end; end
end
