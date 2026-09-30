function report = test_i10p4_m32_restart_reliability()
%TEST_I10P4_M32_RESTART_RELIABILITY  Structural regression for I.10.4.

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.4 - restart reliability regression suite\n');
    fprintf('============================================================\n');

    tests={@test_nested_restart_smoke};
    names={'Nested candidate-vs-extended restart comparison is reproducible'};
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
        error('test_i10p4_m32_restart_reliability:Failure','Commit I.10.4 regression suite failed.');
    end
end

function test_nested_restart_smoke()
    root=tempname; mkdir(root); cleanup=onCleanup(@()safe_rmdir(root)); %#ok<NASGU>
    d=sim_i10p4_m32_restart_reliability( ...
        'M',4,'SNRdB',20,'TurbulenceVec',0.1, ...
        'CandidateStarts',2,'ExtendedStarts',3,'nReplicates',1,'saMaxIter',1, ...
        'UseParallelCases',false,'WriteCSV',false,'ResultsRoot',root,'RunLabel','regression');
    T=d.results;
    assert(height(T)==1,'Expected exactly one smoke row.');
    assert(T.BestExtendedAMI>=T.BestCandidateAMI-1e-14, ...
        'Nested extended winner must not be worse than candidate-prefix winner.');
    assert(T.ExtraValidatedGain>=-1e-14,'Nested extra gain must be nonnegative.');
    assert(T.BestCandidateStart<=2 && T.BestExtendedStart<=3,'Unexpected restart index.');
    assert(abs(T.SelectedDt-0.01)<1e-14 && ~T.IsRefined, ...
        'M4/SNR20 smoke point should remain a base-dt case.');
    assert(T.RuleVersion=="I10.3-v1",'Unexpected dt selector version.');
    assert(d.summary.nRows==1 && d.summary.candidateStarts==2 && d.summary.extendedStarts==3);
    assert(isfinite(d.summary.meanExtendedOverCandidateRuntimeRatio) && ...
        d.summary.meanExtendedOverCandidateRuntimeRatio>0,'Invalid runtime ratio.');
end

function safe_rmdir(p)
    if exist(p,'dir'),try,rmdir(p,'s');catch,end,end
end
