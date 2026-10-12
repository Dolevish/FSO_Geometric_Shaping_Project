function report = test_k4p2_min_gap(varargin)
%TEST_K4P2_MIN_GAP  5 preflight tests; NEVER runs SA or AMI evaluation.
%
%   test_k4p2_min_gap();
%   test_k4p2_min_gap('K4CombinedPath','/path/K4_combined_results.mat');
%
% Test 5 clones original .005 source rows to .001 in temporary MEMORY
% to exercise the actual merge/plot, explicitly as SYNTHETIC fixture.
% It NEVER overwrites or alters any real previous results.

    p=inputParser;p.FunctionName=mfilename;
    addParameter(p,'K4CombinedPath','',@is_text);
    parse(p,varargin{:});o=p.Results;

    fprintf('\n============================================================\n');
    fprintf('COMMIT K4.2 - MINIMUM SPACING FROZEN PREFLIGHT\n');
    fprintf('============================================================\n');
    tests={@check_config,@check_paired_starts,@check_dry_run, ...
        @check_geometry,@()synthetic_merge(o.K4CombinedPath)};
    labels={'Frozen scientific config and 32-task budget', ...
        'Deterministic seeds/start identities compatible with K/K4', ...
        'K4.2 2 points / 4 cases / 32 unique restarts dry-run', ...
        'Feasible d_min=0.001 geometry at M16 and M32', ...
        'SYNTHETIC merge from original .005 source and no-errorbar plot'};
    status=false(5,1);info=strings(5,1);
    for i=1:5
        ticId=tic;
        try
            tests{i}();
            status(i)=true;
            fprintf('[%d/5] %s ... PASS (%.2fs)\n',i,labels{i},toc(ticId));
        catch ME
            info(i)=string(ME.message);
            fprintf('[%d/5] %s ... FAIL (%s)\n',i,labels{i},ME.message);
        end
    end
    report=struct('passed',status,'messages',info,'nPassed',nnz(status), ...
        'nTotal',5);
    fprintf('K4.2 PREFLIGHT: %d/5 PASS\n',nnz(status));
    if ~all(status)
        error('K4p2:Test','Preflight failed; do not start production.');
    end
end

function check_config()
    A=k_min_gap_config();K=k4p2_min_gap_config();
    assert(strcmp(A.version,'K-v1') && strcmp(K.version,'K4p2-v1'));
    assert(isequal(K.MinGapVec,.001) && ...
        isequal(K.MinGapByM,{.001,.001}));
    assert(isequal(K.MVec,[16 32]) && K.SNRdB==20 && K.SigmaR2==.2);
    assert(K.nPhysicalSensitivityPoints==2 && ...
        K.nReplicateCases==4 && K.nRestartTasks==32);
    assert(K.totalSAIterations==256000 && ...
        K.nominalFastObjectiveEvaluations==256032);
    fields={'MVec','SNRdB','SigmaR2','Replicates','NStarts','MaxIter', ...
        'P_avg','pinZero','baseSeed','FastFadingMethod', ...
        'LoghDtPolicy','loghDt','loghDtRefined','SelectedDtByM', ...
        'loghSpanSigma','loghYBlockSize','YGridMode','YGridScale', ...
        'historyEvery','T0','Tf','BaseStd0','ItersPerTemp', ...
        'scheduler','taskGranularity','longestFirst','noNestedParallelism'};
    for i=1:numel(fields)
        f=fields{i};
        assert(isequaln(A.(f),K.(f)),'Changed policy field: %s',f);
    end
end

function check_paired_starts()
    A=k_min_gap_config(); E=k4_min_gap_extension_config();
    K=k4p2_min_gap_config();
    [CA,~]=k_build_task_plan(A);
    [CE,~]=k_build_task_plan(E);
    [CK,TK]=k_build_task_plan(K);
    assert(numel(unique(TK.IdentityKey))==32);
    for m=[16 32]
        for rep=[1 2]
            a=unique(CA.Seed(CA.M==m & CA.Replicate==rep));
            b=unique(CE.Seed(CE.M==m & CE.Replicate==rep));
            c=unique(CK.Seed(CK.M==m & CK.Replicate==rep));
            assert(numel(a)==1 && numel(b)==1 && numel(c)==1 ...
                && a==b && b==c);
        end
    end
end

function check_dry_run()
    txt=evalc("D=sim_k4p2_min_gap('DryRun',true);");
    assert(D.dryRun && strcmp(D.version,'K4p2-v1'));
    assert(height(D.casePlan)==4 && height(D.taskPlan)==32);
    assert(all(abs(D.casePlan.MinGap-.001)<1e-12));
    assert(contains(txt,'restart tasks=32'));
end

function check_geometry()
    K=k4p2_min_gap_config();
    for m=K.MVec
        x=define_constellation(m,K.P_avg,true,"mean");
        assert(abs(mean(x)-K.P_avg)<1e-10 && abs(x(1))<1e-12);
        assert(all(diff(x)>=.001-1e-12));
        assert(.001<2*K.P_avg/(m-1));
    end
end

function synthetic_merge(given)
    % Read only REAL original K4 data, but create a FAKE added result in
    % memory for testing code. Must never publish synthetic AMI.
    if strlength(string(given))>0
        source=char(string(given));
    else
        base=fso_result_utils.revision2_results_root();
        D=dir(fullfile(base,'commit_k_min_gap_extension', ...
            'K_*','combined','K4_combined_results.mat'));
        if numel(D)~=1
            error('K4p2:TestSource', ...
                ['Expected one K4 combined MAT under results/ieee_revision2. ' ...
                 'Found %d; provide K4CombinedPath explicitly.'],numel(D));
        end
        source=fullfile(D(1).folder,D(1).name);
    end
    S=load(source,'result');
    if ~isfield(S,'result'),error('K4p2:TestSource','No result in source MAT.');end
    original=S.result;
    K=k4p2_min_gap_config();
    dummy=struct('version','K4p2-v1','config',K, ...
        'baselines',original.baselines, ...
        'outputDirectory',fso_result_utils.revision2_results_root());
    old=.005;
    rt=original.restartSummary(abs(original.restartSummary.MinGap-old)<1e-12,:);
    rp=original.replicateSummary(abs(original.replicateSummary.MinGap-old)<1e-12,:);
    ph=original.physicalSummary(abs(original.physicalSummary.MinGap-old)<1e-12,:);
    assert(height(rt)==32 && height(rp)==4 && height(ph)==2);
    rt.MinGap(:)=.001;rp.MinGap(:)=.001;ph.MinGap(:)=.001;
    rt.IdentityKey=replace(rt.IdentityKey,'gap0p005','gap0p001');
    ph.GainRetentionVsGap0(:)=NaN;
    dummy.restartSummary=rt;
    dummy.replicateSummary=rp;
    dummy.physicalSummary=ph;

    root=fso_result_utils.revision2_results_root();
    temp=tempname(root);
    mkdir(temp);
    cleanup=onCleanup(@() remove_temp(temp)); %#ok<NASGU>
    C=merge_k4p2_min_gap(original,dummy,'OutputDirectory',temp);
    assert(height(C.restartSummary)==256 && ...
        height(C.replicateSummary)==32 && ...
        height(C.physicalSummary)==16);
    assert(~any(abs(C.physicalSummary.MinGap)<1e-12));
    assert(nnz(abs(C.physicalSummary.MinGap-.001)<1e-12)==2);
    assert(exist(C.figurePaths.pdf,'file')==2 && ...
        exist(C.figurePaths.png,'file')==2 && ...
        exist(C.figurePaths.fig,'file')==2);
    assert(all(abs(C.physicalSummary.GainRetentionVsGap0001( ...
        abs(C.physicalSummary.MinGap-.001)<1e-12)-1)<1e-12));
    assert(~ismember('GainRetentionVsGap0', ...
        C.physicalSummary.Properties.VariableNames));
    % Cannot call the synthetic output REAL production data.
    fprintf('  Synthetic fixture only (NOT new optimizer results).\n');
    close all;
end
function remove_temp(d)
    if exist(d,'dir')==7,try,rmdir(d,'s');catch,end,end
end
function tf=is_text(x),tf=ischar(x)||(isstring(x)&&isscalar(x));end
