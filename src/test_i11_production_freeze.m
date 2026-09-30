function report = test_i11_production_freeze()
%TEST_I11_PRODUCTION_FREEZE  Regression checks for I.11 freeze/resume logic.

    fprintf('\n============================================================\n');
    fprintf('Commit I.11 - production freeze regression suite\n');
    fprintf('============================================================\n');

    tests={@test_frozen_policy,@test_variable_restart_scheduler, ...
        @test_checkpoint_resume,@test_resume_mismatch_rejected};
    names={'Frozen scientific profile and restart map are exact', ...
        'Flat scheduler executes variable 6/8 restart counts by M', ...
        'Completed restart checkpoints are reused exactly on resume', ...
        'Resume refuses a mismatched scientific signature'};

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
        error('test_i11_production_freeze:Failure','I.11 production-freeze regression suite failed.');
    end
end

function test_frozen_policy()
    F=i11_production_config();
    assert(strcmp(F.freezeVersion,'I11-v1'));
    assert(isequal(F.MVec,[4 8 16 32]));
    assert(isequal(F.SNRVec,5:5:30));
    assert(max(abs(F.TurbulenceVec-[0 .1 .2 .3]))<1e-14);
    assert(F.saMaxIter==4000 && F.nReplicates==2);
    assert(abs(F.T0-0.4)<1e-14 && abs(F.Tf-1e-3)<1e-14 && abs(F.BaseStd0-0.15)<1e-14);
    assert(F.ItersPerTemp==50 && abs(F.minGap-0.01)<1e-14 && F.pinZero);
    assert(strcmp(F.dtRuleVersion,'I10.3-v1') && strcmp(F.runtimeTuningVersion,'I10.2.2-v1'));
    got=zeros(size(F.MVec));
    for k=1:numel(F.MVec)
        [got(k),m]=i11_restart_policy(F.MVec(k));
        assert(strcmp(m.ruleVersion,'I11-v1'));
    end
    assert(isequal(got,[6 6 8 8]));
    assert(isequal(F.RestartsByM,got));
    assert(F.nOptimizationCases==192);
    assert(F.totalRestartTasks==1344);
    assert(~F.globalOptimumClaim);
end

function test_variable_restart_scheduler()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    b=sim_revision_AMI_vs_SNR_checkpointed_flat( ...
        'MVec',[4 16],'SNRVec',10,'TurbulenceVec',0, ...
        'saMaxIter',1,'RestartPolicy','i11','nReplicates',1, ...
        'UseParallel',false,'Checkpoint',false,'WriteCSV',false, ...
        'ResultsRoot',root,'RunLabel','I11_variable_restart_regression');
    assert(isequal(b.sa.startsByM,[6 8]));
    assert(b.parallel.nRestartTasks==14);
    T=b.summary.caseTable;
    assert(T.NStarts(T.M==4)==6);
    assert(T.NStarts(T.M==16)==8);
end

function test_checkpoint_resume()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    common={'MVec',4,'SNRVec',20,'TurbulenceVec',0.1,'saMaxIter',1, ...
        'RestartPolicy','fixed','saNStarts',2,'nReplicates',1, ...
        'UseParallel',false,'Checkpoint',true,'WriteCSV',false, ...
        'ResultsRoot',root,'RunLabel','I11_resume_regression'};

    b1=sim_revision_AMI_vs_SNR_checkpointed_flat(common{:});
    assert(b1.checkpoint.nExecuted==2 && b1.checkpoint.nReused==0);
    assert(~b1.checkpoint.baselineReused);

    b2=sim_revision_AMI_vs_SNR_checkpointed_flat(common{:}, ...
        'ResumeDirectory',b1.run.outputDirectory);
    assert(b2.checkpoint.nExecuted==0 && b2.checkpoint.nReused==2, ...
        'All completed restart tasks must be reused on resume.');
    assert(b2.checkpoint.baselineReused,'Baseline cache should be reused on resume.');
    assert(abs(b1.summary.caseTable.AMIValidated(1)-b2.summary.caseTable.AMIValidated(1))<1e-14);
    assert(b1.summary.caseTable.Seed(1)==b2.summary.caseTable.Seed(1));
end

function test_resume_mismatch_rejected()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    common={'MVec',4,'SNRVec',20,'TurbulenceVec',0.1,'saMaxIter',1, ...
        'RestartPolicy','fixed','saNStarts',1,'nReplicates',1, ...
        'UseParallel',false,'Checkpoint',true,'WriteCSV',false, ...
        'ResultsRoot',root,'RunLabel','I11_mismatch_regression'};
    b=sim_revision_AMI_vs_SNR_checkpointed_flat(common{:});

    threw=false;
    try
        sim_revision_AMI_vs_SNR_checkpointed_flat( ...
            'MVec',4,'SNRVec',20,'TurbulenceVec',0.1,'saMaxIter',2, ...
            'RestartPolicy','fixed','saNStarts',1,'nReplicates',1, ...
            'UseParallel',false,'Checkpoint',true,'WriteCSV',false, ...
            'ResultsRoot',root,'RunLabel','I11_mismatch_regression', ...
            'ResumeDirectory',b.run.outputDirectory);
    catch ME
        threw=strcmp(ME.identifier,'sim_revision_AMI_vs_SNR_checkpointed_flat:ResumeConfigMismatch');
    end
    assert(threw,'A changed scientific configuration must be rejected during resume.');
end

function safe_rmdir(path)
    if exist(path,'dir')==7
        try,rmdir(path,'s');catch,end
    end
end
