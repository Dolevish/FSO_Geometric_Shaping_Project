function report = test_i9_dynamic_logh_dt()
%TEST_I9_DYNAMIC_LOGH_DT  Regression checks for dynamic log-h dt policy.
%
% Historical Commit I.9 introduced the deterministic per-case selector.
% Commit I.10.3 conservatively expands the refined region by adding the
% (M,SNR,sigma_X^2)=(16,30,0.2) boundary case identified by the I.10 pilot.

    fprintf('\n============================================================\n');
    fprintf('Commit I.9 / I.10.3 - dynamic log-h dt regression suite\n');
    fprintf('============================================================\n');

    tests={@test_selector_exact_current_grid,@test_fixed_override, ...
        @test_production_driver_new_boundary_case,@test_policy_sweep_smoke, ...
        @test_canonical_builder_unchanged};
    names={'Selector maps current production grid to exactly five refined cases', ...
        'Fixed override is explicit and deterministic', ...
        'Production driver stores revised M16 boundary dt', ...
        'Dynamic-dt no-SA policy sweep smoke', ...
        'Canonical build_fso_config defaults remain unchanged'};

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
    if ~all(passed), error('test_i9_dynamic_logh_dt:Failure','Dynamic log-h dt regression suite failed.'); end
end

function test_selector_exact_current_grid()
    Mv=[4 8 16 32]; Sv=5:5:30; Gv=[0 0.1 0.2 0.3]; refined=[];
    for M=Mv
        for s=Sv
            for g=Gv
                [dt,m]=select_logh_dt(M,s,g);
                assert(abs(dt-(0.005*m.isRefined+0.01*(~m.isRefined)))<1e-14);
                assert(strcmp(m.ruleVersion,'I10.3-v1'),'Unexpected dt-policy rule version.');
                [dt2,m2]=select_logh_dt(M,s,g);
                assert(dt==dt2 && m.isRefined==m2.isRefined,'Selector must be deterministic.');
                if m.isRefined, refined(end+1,:)=[M s g]; end %#ok<AGROW>
            end
        end
    end
    expected=[16 30 0.2;16 30 0.3;32 30 0.1;32 30 0.2;32 30 0.3];
    refined=sortrows(refined,[1 2 3]); expected=sortrows(expected,[1 2 3]);
    assert(isequal(size(refined),size(expected)) && max(abs(refined(:)-expected(:)))<1e-12, ...
        'Current production grid must contain exactly the five I.10.3 refined cases.');
end

function test_fixed_override()
    [a,ma]=select_logh_dt(32,30,0.3,'Policy','fixed','FixedDt',0.01);
    [b,mb]=select_logh_dt(4,5,0,'Policy','fixed','FixedDt',0.005);
    assert(abs(a-0.01)<1e-14 && ~ma.isRefined && strcmp(ma.policy,'fixed'));
    assert(abs(b-0.005)<1e-14 && ~mb.isRefined && strcmp(mb.policy,'fixed'));
    assert(strcmp(ma.ruleVersion,'I10.3-v1') && strcmp(mb.ruleVersion,'I10.3-v1'));
end

function test_production_driver_new_boundary_case()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    b=sim_revision_AMI_vs_SNR('MVec',16,'SNRVec',30,'TurbulenceVec',0.2, ...
        'saMaxIter',1,'saNStarts',1,'nReplicates',1,'UseParallelCases',false, ...
        'MakePlots',false,'SaveFigures',false,'ResultsRoot',root,'RunLabel','I10p3_regression');
    assert(strcmp(b.method.loghDtPolicy,'case-adaptive'));
    assert(abs(b.summary.loghDtSelected(1)-0.005)<1e-14);
    assert(b.summary.loghDtIsRefined(1));
    p=b.caseFiles{1,1,1,1}; s=load(p,'c');
    assert(abs(s.c.integration.loghDtSelected-0.005)<1e-14);
    assert(s.c.integration.loghDtIsRefined);
    assert(strcmp(s.c.integration.loghDtRuleVersion,'I10.3-v1'));
    assert(abs(s.c.config.loghDt-0.005)<1e-14);
end

function test_policy_sweep_smoke()
    root=tempname; mkdir(root); cleanup=onCleanup(@() safe_rmdir(root)); %#ok<NASGU>
    d=diagnose_production_dt_policy_accuracy('MVec',[4 16],'SNRVec',[20 30], ...
        'TurbulenceVec',0.2,'GeometryModes',{'uniform'},'UseParallelCases',false, ...
        'WriteCSV',false,'ResultsRoot',root,'RunLabel','I10p3_diag_regression');
    T=d.results;
    assert(height(T)==4,'Expected one geometry for four physical points.');
    hard=T.M==16 & T.SNRdB==30 & abs(T.SigmaX2-0.2)<1e-12;
    assert(sum(hard)==1 && abs(T.SelectedDt(hard)-0.005)<1e-14 && T.IsRefined(hard));
    assert(all(abs(T.SelectedDt(~hard)-0.01)<1e-14));
    assert(d.summary.nRefinedPhysicalCases==1);
end

function test_canonical_builder_unchanged()
    cfg=build_fso_config(4,1,20,0.1,'minGap',0.01,'saMaxIter',1,'saNStarts',1);
    assert(strcmp(cfg.yGridMode,'fixed-bound'));
    assert(abs(cfg.loghDt-0.01)<1e-14, ...
        'Canonical build_fso_config default must remain historical dt=0.01.');
end

function safe_rmdir(p)
    if exist(p,'dir'), try, rmdir(p,'s'); catch, end; end
end
