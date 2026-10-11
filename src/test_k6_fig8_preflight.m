function report = test_k6_fig8_preflight(varargin)
%TEST_K6_FIG8_PREFLIGHT  Four science/reproducibility checks + plot smoke.
%
% No SA iterations, no fast objective evaluations, and no independent
% AMI validations are performed by this test. Synthetic values are
% used ONLY for exercising plot output, never as scientific results.
%
%   report=test_k6_fig8_preflight();
%   report=test_k6_fig8_preflight( ...
%       'ProductionBundlePath','/path/i11p1ProductionBundle.mat');

    fprintf('\n============================================================\n');
    fprintf('COMMIT K6 - Fig. 8 frozen I.11.1 preflight\n');
    fprintf('============================================================\n');

    tests={@check_frozen_scientific_policy, ...
        @()check_real_production_plan(varargin{:}), ...
        @check_independent_threefry_initialization, ...
        @check_result_directory, ...
        @check_synthetic_plot_export};
    names={'Original I.11.1 channel, constraints and SA parameters', ...
        'Real production source and 20-restart dry run', ...
        '20 unique Threefry trajectories and exact initial-point policy', ...
        'Results path remains under ieee_revision2', ...
        'Fig. 8 bar-chart PNG/PDF/FIG export (synthetic only)'};

    passed=false(numel(tests),1);
    messages=strings(numel(tests),1);
    durations=zeros(numel(tests),1);
    for k=1:numel(tests)
        ticId=tic;
        try
            tests{k}();
            passed(k)=true;
            messages(k)="PASS";
            fprintf('[%d/%d] %s ... PASS (%.2f s)\n', ...
                k,numel(tests),names{k},toc(ticId));
        catch ME
            messages(k)=string(ME.message);
            fprintf('[%d/%d] %s ... FAIL (%s)\n', ...
                k,numel(tests),names{k},ME.message);
        end
        durations(k)=toc(ticId);
    end
    report=struct('passed',passed,'messages',messages, ...
        'elapsedSeconds',durations, ...
        'nPassed',nnz(passed),'nTotal',numel(tests));
    fprintf('\nK6 PREFLIGHT %d/%d PASS\n',report.nPassed,report.nTotal);
    if ~all(passed)
        error('test_k6_fig8_preflight:Failed', ...
            'Do not start K6 until all preflight tests pass.');
    end
end


function check_frozen_scientific_policy()
    F=i11p1_production_config();
    assert(strcmp(F.freezeVersion,'I11.1-v1'));
    assert(F.P_avg==1 && F.minGap==0.01 && F.pinZero);
    assert(isequal(F.MVec,[4 8 16 32]));
    assert(isequal(F.RestartsByM,[6 6 8 8]));
    assert(isequal(F.MaxIterByM,[4000 4000 4000 8000]));
    assert(i11p1_iteration_policy(8)==4000);
    assert(isequal(F.TurbulenceVec,[0 0.1 0.2 0.3]));
    assert(isequal(F.SNRVec,5:5:30));
    assert(strcmp(F.FastFadingMethod,'logh'));
    assert(strcmp(F.YGridMode,'candidate-adaptive'));
    assert(F.loghDt==0.01 && F.loghDtRefined==0.005);
    assert(F.YGridScale==1.1 && F.loghSpanSigma==8);
    assert(F.loghYBlockSize==512);
    assert(F.T0==0.4 && F.Tf==0.001);
    assert(F.BaseStd0==0.15 && F.ItersPerTemp==50);
    assert(F.historyEvery==50 && F.baseSeed==20260928);
end


function check_real_production_plan(varargin)
    d=sim_k6_fig8_repeatability('DryRun',true,varargin{:});
    assert(d.dryRun && strcmp(d.version,'K6-FIG8-I11p1-v1'));
    assert(d.nRestarts==20 && d.totalSAIterations==80000);
    assert(height(d.taskPlan)==20 && ...
        numel(unique(d.taskPlan.IdentityKey))==20);
    assert(isequal(d.taskPlan.Restart(:),(1:20).'));
    assert(isequal(d.taskPlan.Substream(:),(1:20).'));
    assert(numel(unique(d.taskPlan.MasterSeed))==1);
    assert(height(d.sourceCases)==2 && ...
        isequal(d.sourceCases.Replicate(:).',[1 2]));
    assert(~any(d.taskPlan.MasterSeed(1)==d.sourceCases.Seed));
    C=d.config;
    assert(C.M==8 && C.P_avg==1 && C.SNR_dB==20);
    assert(abs(C.sigma_X_sq-.1)<1e-12);
    assert(C.SA.maxIter==4000 && C.SA.nStarts==20);
    assert(C.SA.minGap==.01 && C.SA.pinZero);
    assert(C.SA.historyEvery==50 && C.SA.logEvery==0);
    assert(strcmp(C.fastFadingMethod,'logh'));
    assert(strcmp(C.yGridMode,'candidate-adaptive'));
    assert(abs(C.loghDt-.01)<1e-12);
    assert(isequal(C.px,ones(8,1)/8));
    assert(C.xMaxBound==fso_result_utils.feasible_xmax_bound(8,1,.01,true));
end


function check_independent_threefry_initialization()
    F=i11p1_production_config();
    seed=uint32(fso_result_utils.case_seed(20261011,8,20,.1));
    b=define_constellation(8,F.P_avg,true,"mean");
    cfg=build_fso_config(8,1,20,.1, ...
        'FastFadingMethod',F.FastFadingMethod, ...
        'loghDt',.01,'loghSpanSigma',F.loghSpanSigma, ...
        'loghYBlockSize',F.loghYBlockSize, ...
        'YGridMode',F.YGridMode,'YGridScale',F.YGridScale, ...
        'xMaxBound',fso_result_utils.feasible_xmax_bound(8,1,.01,true), ...
        'saMaxIter',4000,'saNStarts',20, ...
        'minGap',F.minGap,'pinZero',F.pinZero, ...
        'seedInit',double(seed),'logEvery',0,'historyEvery',50);
    X=i13_shared_start_points(cfg,b,seed,20);
    assert(size(X,2)==20);
    assert(isequal(X(:,1),b(:)));
    assert(isequal(X,i13_shared_start_points(cfg,b,seed,20)));
    for k=2:20
        s=RandStream('Threefry','Seed',double(seed));
        s.Substream=k;
        expected=b(:)+.2*randn(s,8,1);
        assert(isequal(X(:,k),expected), ...
            'K6 initializations must reproduce the frozen SA restart rule.');
    end
    assert(size(unique(X(:,2:end).','rows'),1)==19);
    assert(exist('sa_restart_task_checkpointed','file')==2);
    assert(exist('sa_restart_task','file')==2);
end


function check_result_directory()
    root=fso_result_utils.revision2_results_root();
    [~,name]=fileparts(root);
    assert(strcmp(name,'ieee_revision2'));
    assert(~strcmp(root,fso_result_utils.default_results_root()));
end


function check_synthetic_plot_export()
    root=fso_result_utils.revision2_results_root();
    if exist(root,'dir')~=7,mkdir(root);end
    tmp=tempname(root);
    mkdir(tmp);
    cleanup=onCleanup(@() cleanup_dir(tmp)); %#ok<NASGU>

    % Exercise narrow-spread inset with SYNTHETIC mock results.
    y=1.695+0.0003*(1:20).' + ...
        0.0005*sin((1:20).'*1.7);
    T=table((1:20).',y, ...
        'VariableNames',{'Restart','ValidatedAMI'});
    S=summary_for_test(y);
    paths=plot_k6_fig8_repeatability(T,S, ...
        'OutputDirectory',tmp,'CloseFigure',true);
    assert(exist(paths.pdf,'file')==2);
    assert(exist(paths.png,'file')==2);
    assert(exist(paths.fig,'file')==2);

    shouldFail=false;
    S.MeanAMI=S.MeanAMI+0.01;
    try
        plot_k6_fig8_repeatability(T,S, ...
            'OutputDirectory',tmp,'CloseFigure',true);
    catch
        shouldFail=true;
    end
    assert(shouldFail,'Figure must reject inconsistent summary statistics.');
end


function S=summary_for_test(y)
    S=struct();
    S.MeanAMI=mean(y);
    S.SampleVarianceAMI=var(y,0);
    S.SampleStdAMI=std(y,0);
    S.MinAMI=min(y);
    S.MaxAMI=max(y);
    S.RangeAMI=max(y)-min(y);
end


function cleanup_dir(path)
    if exist(path,'dir')==7,try,rmdir(path,'s');catch,end,end
end
