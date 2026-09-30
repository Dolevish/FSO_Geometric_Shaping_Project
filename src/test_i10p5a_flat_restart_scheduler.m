function report = test_i10p5a_flat_restart_scheduler()
%TEST_I10P5A_FLAT_RESTART_SCHEDULER  Regression suite for flat restart tasks.

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.5A - flat restart-task scheduler regression suite\n');
    fprintf('============================================================\n');

    tests={@test_atomic_restart_matches_sa_multistart, ...
           @test_flat_driver_matches_legacy_serial, ...
           @test_parallel_flat_smoke_if_available};
    names={'Atomic restart reconstruction matches legacy sa_multistart', ...
           'Flat serial driver matches legacy case driver', ...
           'Parallel flat queue smoke (when Parallel Computing Toolbox is available)'};

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
        error('test_i10p5a_flat_restart_scheduler:Failure','I.10.5A regression suite failed.');
    end
end


function test_atomic_restart_matches_sa_multistart()
    M=4; Pavg=1; snr=20; sig=0.1; seed=239232940;
    cfg=build_fso_config(M,Pavg,snr,sig, ...
        'FastFadingMethod','logh','loghDt',0.01,'loghSpanSigma',8, ...
        'YGridMode','candidate-adaptive','YGridScale',1.1, ...
        'saMaxIter',8,'saNStarts',3,'saUseParallel',false, ...
        'minGap',0.01,'pinZero',true,'seedInit',seed,'logEvery',0,'historyEvery',2);
    cfg.SA.T0=0.4; cfg.SA.Tf=1e-3; cfg.SA.baseStd0=0.15; cfg.SA.itersPerTemp=4;
    cfg.SA.nBlocks=ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate=exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);

    x=define_constellation(M,Pavg,true,"mean");
    [oldOut,oldStarts]=sa_multistart(cfg,x,'[I10p5A-legacy]');

    newStarts=cell(3,1);
    for r=1:3
        newStarts{r}=sa_restart_task(cfg,x,seed,r,'[I10p5A-flat]');
    end
    [newOut,newStarts]=sa_finalize_restart_results(vertcat(newStarts{:}),'[I10p5A-flat]');

    for r=1:3
        assert(max(abs(oldStarts(r).x0_raw-newStarts(r).x0_raw))<1e-14,'x0_raw mismatch at restart %d.',r);
        assert(max(abs(oldStarts(r).x0_projected-newStarts(r).x0_projected))<1e-14,'x0_projected mismatch.');
        assert(abs(oldStarts(r).bestMIFast-newStarts(r).bestMIFast)<5e-13,'Fast AMI mismatch.');
        assert(abs(oldStarts(r).bestMIValidated-newStarts(r).bestMIValidated)<5e-13,'Validated AMI mismatch.');
        assert(max(abs(oldStarts(r).x_best-newStarts(r).x_best))<1e-13,'Winner geometry mismatch.');
        assert(isequal(oldStarts(r).history.iter,newStarts(r).history.iter),'History iteration mismatch.');
        assert(max(abs(oldStarts(r).history.bestMI-newStarts(r).history.bestMI))<5e-13,'History AMI mismatch.');
    end
    assert(oldOut.bestStartValidated==newOut.bestStartValidated,'Validated winner mismatch.');
    assert(oldOut.bestStartFast==newOut.bestStartFast,'Fast winner mismatch.');
    assert(abs(oldOut.bestMIValidated-newOut.bestMIValidated)<5e-13,'Final validated AMI mismatch.');
end


function test_flat_driver_matches_legacy_serial()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    common={'MVec',4,'SNRVec',20,'TurbulenceVec',0.1,'P_avg',1, ...
        'minGap',0.01,'pinZero',true,'FastFadingMethod','logh', ...
        'LoghDtPolicy','case-adaptive','loghDt',0.01,'loghDtRefined',0.005, ...
        'loghSpanSigma',8,'YGridMode','candidate-adaptive','YGridScale',1.1, ...
        'saMaxIter',6,'saNStarts',3,'nReplicates',1, ...
        'T0',0.4,'Tf',1e-3,'BaseStd0',0.15,'ItersPerTemp',3, ...
        'baseSeed',20260928,'historyEvery',2};

    legacy=sim_revision_AMI_vs_SNR(common{:}, ...
        'UseParallelCases',false,'MakePlots',false,'SaveFigures',false, ...
        'ResultsRoot',fullfile(root,'legacy'),'RunLabel','legacy');

    flat=sim_revision_AMI_vs_SNR_flat_restarts(common{:}, ...
        'UseParallel',false,'WriteCSV',false, ...
        'ResultsRoot',fullfile(root,'flat'),'RunLabel','flat');

    a=load(legacy.caseFiles{1,1,1,1},'c'); b=load(flat.caseFiles{1,1,1,1},'c');
    ca=a.c; cb=b.c;
    assert(abs(ca.integration.loghDtSelected-cb.integration.loghDtSelected)<1e-14);
    assert(ca.winner.startIndex==cb.winner.startIndex,'Winner restart mismatch.');
    assert(abs(ca.winner.amiValidated-cb.winner.amiValidated)<5e-13,'Validated winner mismatch.');
    assert(abs(ca.winner.amiFast-cb.winner.amiFast)<5e-13,'Fast winner mismatch.');
    assert(max(abs(ca.winner.x-cb.winner.x))<1e-13,'Winner constellation mismatch.');
    assert(numel(ca.starts)==numel(cb.starts));
    for r=1:numel(ca.starts)
        assert(abs(ca.starts(r).bestMIFast-cb.starts(r).bestMIFast)<5e-13);
        assert(abs(ca.starts(r).bestMIValidated-cb.starts(r).bestMIValidated)<5e-13);
        assert(max(abs(ca.starts(r).x_best-cb.starts(r).x_best))<1e-13);
    end
    assert(strcmp(flat.parallel.granularity,'restart-task') && flat.parallel.noNestedParallelism);
end


function test_parallel_flat_smoke_if_available()
    hasPCT=false;
    try, hasPCT=license('test','Distrib_Computing_Toolbox')~=0; catch, end
    if ~hasPCT
        fprintf('[PCT unavailable: structural skip] ');
        return;
    end

    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    b=sim_revision_AMI_vs_SNR_flat_restarts( ...
        'MVec',4,'SNRVec',20,'TurbulenceVec',[0.1 0.2], ...
        'saMaxIter',1,'saNStarts',3,'nReplicates',1, ...
        'ItersPerTemp',1,'historyEvery',1, ...
        'UseParallel',true,'LongestFirst',true,'WriteCSV',false, ...
        'ResultsRoot',root,'RunLabel','parallel_smoke');

    assert(b.parallel.useParallel,'Expected parallel execution when PCT is available.');
    assert(b.parallel.nRestartTasks==6,'Expected 2 cases x 3 restart tasks.');
    assert(b.parallel.noNestedParallelism);
    assert(all(isfinite(b.summary.caseTable.AMIValidated)));
end


function safe_rmdir(p)
    if exist(p,'dir'), try, rmdir(p,'s'); catch, end; end
end
