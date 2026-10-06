function report = test_k_min_gap_sensitivity()
%TEST_K_MIN_GAP_SENSITIVITY  Structural preflight for Commit K.
%
% This test does not run SA. It verifies the frozen sensitivity grid,
% paired-seed policy, evaluator resolution, task count, feasibility and the
% new revision-2 results-root convention.

    fprintf('\n============================================================\n');
    fprintf('COMMIT K - minimum-spacing sensitivity preflight\n');
    fprintf('============================================================\n');

    tests={@test_config,@test_task_plan,@test_paired_seeds,@test_case_cfgs,@test_root};
    names={ ...
        'Frozen grid and compute budget', ...
        '20 replicate cases / 160 unique restart tasks', ...
        'Same master seed across d_min for each M/replicate', ...
        'Production evaluator and all d_min configurations are feasible', ...
        'Results root is results/ieee_revision2'};

    passed=false(numel(tests),1); elapsed=nan(numel(tests),1); messages=strings(numel(tests),1);
    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k});
        t=tic;
        try
            tests{k}();
            passed(k)=true; messages(k)="PASS"; elapsed(k)=toc(t);
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
        error('test_k_min_gap_sensitivity:Failure','Commit-K preflight failed.');
    end
end


function test_config()
    K=k_min_gap_config();
    assert(isequal(K.MVec,[16 32]));
    assert(K.SNRdB==20);
    assert(abs(K.SigmaR2-0.2)<eps);
    assert(isequal(K.MinGapVec,[0 0.005 0.01 0.02 0.05]));
    assert(isequal(K.Replicates,[1 2]));
    assert(K.NStarts==8);
    assert(K.MaxIter==8000);
    assert(K.nPhysicalSensitivityPoints==10);
    assert(K.nReplicateCases==20);
    assert(K.nRestartTasks==160);
    assert(K.totalSAIterations==1280000);
    assert(K.nominalFastObjectiveEvaluations==1280160);
    assert(all(abs(K.SelectedDtByM-0.01)<1e-12));
    assert(strcmp(K.FastFadingMethod,'logh'));
    assert(strcmp(K.YGridMode,'candidate-adaptive'));
    assert(abs(K.YGridScale-1.1)<eps);
end


function test_task_plan()
    K=k_min_gap_config();
    [C,Q]=k_build_task_plan(K);
    assert(height(C)==20);
    assert(height(Q)==160);
    assert(numel(unique(Q.IdentityKey))==160);

    for m=K.MVec
        for d=K.MinGapVec
            for rep=K.Replicates
                idx=C.M==m & abs(C.MinGap-d)<1e-12 & C.Replicate==rep;
                assert(sum(idx)==1);
                assert(C.NStarts(idx)==8);
                assert(C.MaxIter(idx)==8000);
                assert(abs(C.SelectedDt(idx)-0.01)<1e-12);
            end
        end
    end

    % Longest-first: M32/d_min=0 tasks should lead.
    first=C(Q.TaskCaseIndex(1),:);
    assert(first.M==32);
    assert(abs(first.MinGap)<1e-12);
end


function test_paired_seeds()
    K=k_min_gap_config();
    [C,~]=k_build_task_plan(K);
    for m=K.MVec
        for rep=K.Replicates
            seeds=C.Seed(C.M==m & C.Replicate==rep);
            assert(numel(seeds)==numel(K.MinGapVec));
            assert(numel(unique(seeds))==1, ...
                'Seed must be paired across all d_min values.');
        end
    end
end


function test_case_cfgs()
    K=k_min_gap_config();
    [C,~]=k_build_task_plan(K);

    for m=K.MVec
        for d=K.MinGapVec
            R=C(find(C.M==m & abs(C.MinGap-d)<1e-12,1),:);
            xmax=fso_result_utils.feasible_xmax_bound(m,K.P_avg,d,K.pinZero);
            cfg=build_fso_config(m,K.P_avg,K.SNRdB,K.SigmaR2, ...
                'FastFadingMethod',K.FastFadingMethod,'loghDt',R.SelectedDt, ...
                'loghSpanSigma',K.loghSpanSigma,'loghYBlockSize',K.loghYBlockSize, ...
                'YGridMode',K.YGridMode,'YGridScale',K.YGridScale,'xMaxBound',xmax, ...
                'saMaxIter',K.MaxIter,'saNStarts',K.NStarts,'saUseParallel',false, ...
                'minGap',d,'pinZero',K.pinZero,'seedInit',R.Seed, ...
                'logEvery',0,'historyEvery',K.historyEvery);

            assert(cfg.SA.maxIter==8000);
            assert(cfg.SA.nStarts==8);
            assert(abs(cfg.SA.minGap-d)<1e-12);
            assert(abs(cfg.loghDt-0.01)<1e-12);
            assert(strcmp(cfg.fastFadingMethod,'logh'));
            assert(strcmp(cfg.yGridMode,'candidate-adaptive'));

            x=define_constellation(m,K.P_avg,true,"mean");
            AMI_functions.assert_constellation_feasible(x,cfg,1e-10);
        end
    end
end


function test_root()
    root=fso_result_utils.revision2_results_root();
    suffix=fullfile('results','ieee_revision2');
    assert(endsWith(root,suffix));
end
