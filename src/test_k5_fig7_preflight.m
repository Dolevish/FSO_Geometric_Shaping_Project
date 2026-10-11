function report = test_k5_fig7_preflight(varargin)
%TEST_K5_FIG7_PREFLIGHT  Source/constraint checks and synthetic figure smoke.
%
% NEVER runs SA, Pattern Search or any AMI evaluator.
% Requires existing frozen I.11.1 production and I.13B Zhou MAT files.
%
%   test_k5_fig7_preflight()
%   test_k5_fig7_preflight('ProductionBundlePath',path, ...
%                         'ZhouResultPath',zhouPath)

    fprintf('\n=======================================================\n');
    fprintf('COMMIT K5 - FIG. 7 FROZEN PROTOCOL PREFLIGHT\n');
    fprintf('=======================================================\n');

    tests={@check_k5_config,@()check_real_sources(varargin{:}), ...
        @check_deterministic_start_policy,@check_synthetic_figure};
    labels={'I.13D constraints, budgets and PS options unchanged', ...
        'Real I.11.1/I.13B sources and 12 restart identities', ...
        'Same Threefry seeded raw restart starts as I.13D', ...
        'Four-method figure exports using SYNTHETIC data'};
    passed=false(numel(tests),1);
    messages=strings(numel(tests),1);
    durations=zeros(numel(tests),1);

    for k=1:numel(tests)
        ti=tic;
        try
            tests{k}();
            passed(k)=true;
            messages(k)="PASS";
            fprintf('[%d/%d] %s ... PASS (%.2f s)\n', ...
                k,numel(tests),labels{k},toc(ti));
        catch ME
            messages(k)=string(ME.message);
            fprintf('[%d/%d] %s ... FAIL (%s)\n', ...
                k,numel(tests),labels{k},ME.message);
        end
        durations(k)=toc(ti);
    end
    report=struct('passed',passed,'messages',messages, ...
        'elapsedSeconds',durations,'nPassed',nnz(passed),'nTotal',numel(tests));
    fprintf('\nK5 PREFLIGHT %d/%d PASS\n',report.nPassed,report.nTotal);
    if ~all(passed)
        error('test_k5_fig7_preflight:Failed', ...
            'Do not start K5 until all four preflight tests pass.');
    end
end


function check_k5_config()
    A=i13d_benchmark_config();
    P=i11p1_production_config();
    assert(A.P_avg==1 && P.P_avg==1 && A.minGap==.01 && P.minGap==.01);
    assert(A.pinZero && P.pinZero);
    assert(isequal(A.Replicates,[1 2]));
    assert(A.RestartsByM(A.MVec==8)==6);
    assert(A.MaxIterByM(A.MVec==8)==4000);
    assert(A.MaxEvalsPerStart(A.MVec==8)==4001);
    assert(~A.PatternSearch.UseCompletePoll && ~A.PatternSearch.UseParallel);
    assert(A.PatternSearch.MeshTolerance==1e-6);
    assert(A.PatternSearch.StepTolerance==1e-8);
    assert(A.PatternSearch.FunctionTolerance==1e-8);
    assert(A.PatternSearch.MaxFinalRepairInfNorm==1e-5);
    assert(strcmp(P.FastFadingMethod,'logh'));
    assert(strcmp(P.YGridMode,'candidate-adaptive'));
    assert(P.YGridScale==1.1 && P.loghSpanSigma==8);
    assert(P.loghYBlockSize==512);
end


function check_real_sources(varargin)
    dry=fig7_m8_20db_sigma01('DryRun',true,varargin{:});
    C=dry.casePlan;T=dry.taskPlan;F=dry.config;
    assert(dry.dryRun);
    assert(height(C)==2 && height(T)==12);
    assert(F.nPhysicalPoints==1 && F.nOptimizationCases==2);
    assert(F.nRestartTasks==12 && F.totalMaxFunctionEvaluations==48012);
    assert(all(C.M==8 & C.SNRdB==20 & abs(C.SigmaR2-.1)<1e-12));
    assert(isequal(sort(C.Replicate).',[1 2]));
    assert(all(C.NStarts==6 & C.MaxIterSA==4000 & ...
        C.MaxEvalsPerStart==4001));
    assert(all(abs(C.SelectedDt-.01)<1e-12));
    assert(numel(unique(T.IdentityKey))==12);
    assert(numel(unique(C.Seed))==2);
    assert(dry.PAM>0 && dry.Zhou>0 && dry.SA>0);
end


function check_deterministic_start_policy()
    F=i13d_benchmark_config();
    cfg=build_fso_config(8,F.P_avg,20,.1, ...
        'FastFadingMethod','logh','loghDt',.01, ...
        'loghSpanSigma',8,'loghYBlockSize',512, ...
        'YGridMode','candidate-adaptive','YGridScale',1.1, ...
        'minGap',F.minGap,'pinZero',F.pinZero, ...
        'saMaxIter',4000,'saNStarts',6, ...
        'seedInit',20260928,'logEvery',0,'historyEvery',50);
    x=define_constellation(8,F.P_avg,true,"mean");
    seed=uint32(20260928);
    X=i13_shared_start_points(cfg,x,seed,6);
    Y=i13_shared_start_points(cfg,x,seed,6);
    assert(isequaln(X,Y) && isequal(X(:,1),x(:)));
    assert(size(X,2)==6);
end


function check_synthetic_figure()
    % Explicitly synthetic values; NO conclusion about Fig.7 performance.
    data=table(["Uniform PAM";"Zhou GS";"Pattern Search";"SA"], ...
        [1.0;1.2;1.25;1.3],[0;0;0.003;0.002], ...
        'VariableNames',{'Method','ValidatedAMI','StdAcrossReplicates'});
    root=fso_result_utils.revision2_results_root();
    if exist(root,'dir')~=7,mkdir(root);end
    tmp=tempname(root);
    mkdir(tmp);
    cleanup=onCleanup(@() cleanup_temp(tmp)); %#ok<NASGU>
    paths=plot_fig7_four_methods(data,.1,tmp,'CloseFigure',true);
    assert(exist(paths.pdf,'file')==2 && ...
        exist(paths.png,'file')==2 && exist(paths.fig,'file')==2);
    got=false;
    try
        plot_fig7_four_methods(data,.3,tmp,'CloseFigure',true);
    catch
        got=true;
    end
    assert(got,'Plot must reject sigma_R^2 different from 0.1.');
end


function cleanup_temp(tmp)
    if exist(tmp,'dir')==7,try,rmdir(tmp,'s');catch,end,end
end
