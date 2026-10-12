function report = test_k5p2_fig7_preflight(varargin)
%TEST_K5P2_FIG7_PREFLIGHT  Verify K5.2 before publication export.
%
% NO optimization, no AMI evaluation. Synthetic plot test is NOT results.
%
%   report=test_k5p2_fig7_preflight();
%   report=test_k5p2_fig7_preflight('I13DResultPath','...', ...
%       'ProductionBundlePath','...', 'ZhouResultPath','...', ...
%       'K5ResultPath','...');

    fprintf('\n==========================================================\n');
    fprintf('COMMIT K5.2 - FIG. 7 FROZEN-RESULT REPLAY PREFLIGHT\n');
    fprintf('==========================================================\n');
    tests={@check_frozen_policy, ...
        @()check_original_and_k5_sources(varargin{:}), ...
        @()check_restarts_and_aggregation(varargin{:}), ...
        @check_synthetic_plot_export, ...
        @check_plot_rejects_bad_science};
    captions={'Frozen original K5/I13D budgets and scientific config', ...
        'I.11.1, I.13B, I.13D and K5 source consistency', ...
        '12 validated I13D restarts and two-replicate reduction', ...
        'SYNTHETIC four-method plot exports PDF/PNG/FIG', ...
        'Plot rejects wrong sigma_R^2 and mismatched methods'};
    passed=false(5,1);errors=strings(5,1);
    for k=1:5
        t=tic;
        try
            tests{k}();
            passed(k)=true;
            fprintf('[%d/5] %s ... PASS (%.2fs)\n',k,captions{k},toc(t));
        catch ME
            errors(k)=string(ME.message);
            fprintf('[%d/5] %s ... FAIL: %s\n',k,captions{k},ME.message);
        end
    end
    report=struct('passed',passed,'errors',errors, ...
        'passedCount',nnz(passed),'total',5);
    fprintf('K5.2 PREFLIGHT: %d/5 PASS\n',report.passedCount);
    if ~all(passed)
        error('K5p2:Preflight','Preflight failed; do not use results for publication.');
    end
end

function check_frozen_policy()
    F=i13d_benchmark_config();
    P=i11p1_production_config();
    assert(F.P_avg==1 && P.P_avg==1 && F.minGap==.01 && P.minGap==.01);
    assert(F.pinZero && P.pinZero && F.tieToleranceBits==1e-3);
    assert(isequal(F.Replicates,[1 2]));
    assert(F.RestartsByM(F.MVec==8)==6);
    assert(F.MaxIterByM(F.MVec==8)==4000);
    assert(F.MaxEvalsPerStart(F.MVec==8)==4001);
    assert(~F.PatternSearch.UseCompletePoll && ...
        ~F.PatternSearch.UseParallel && ...
        F.PatternSearch.MeshTolerance==1e-6 && ...
        F.PatternSearch.StepTolerance==1e-8 && ...
        F.PatternSearch.FunctionTolerance==1e-8 && ...
        F.PatternSearch.MaxFinalRepairInfNorm==1e-5);
    assert(strcmp(P.FastFadingMethod,'logh') && ...
        strcmp(P.YGridMode,'candidate-adaptive'));
    assert(P.YGridScale==1.1 && P.loghSpanSigma==8 && ...
        P.loghYBlockSize==512 && P.historyEvery==50);
end

function check_original_and_k5_sources(varargin)
    D=make_k5p2_fig7('DryRun',true,varargin{:});
    assert(D.dryRun && ~D.optimizationRerun && D.M==8 && ...
        D.SNRdB==20 && abs(D.SigmaR2-.3)<1e-12);
    assert(isequal(string(D.fourMethodTable.Method), ...
        ["Uniform PAM";"Zhou GS";"Pattern Search";"SA"]));
    assert(all(D.fourMethodTable.ValidatedAMI>=0) && ...
        all(D.fourMethodTable.ValidatedAMI<=3));
    assert(all(isfinite(D.k5Comparison.K5_sigmaR2_0p1)));
    assert(all(isfinite(D.k5Comparison.K5p2_sigmaR2_0p3)));
    assert(D.optimizerGapByTurbulence.SigmaR2(1)==0.1 && ...
        D.optimizerGapByTurbulence.SigmaR2(2)==0.3);
end

function check_restarts_and_aggregation(varargin)
    D=make_k5p2_fig7('DryRun',true,varargin{:});
    T=D.frozenRestartRows;
    C=D.frozenCaseRows;
    assert(height(T)==12 && height(C)==2);
    assert(numel(unique(T.IdentityKey))==12);
    assert(all(T.EvalBudget==4001 & T.FunctionEvaluations<=4001));
    assert(all(abs(T.FastValidatorGap)<1e-3));
    assert(all(C.SelectedDt==0.01));
    assert(isequal(sort(C.Replicate).',[1 2]));
    for i=1:2
        A=T(T.Replicate==i,:);
        assert(height(A)==6 && isequal(sort(A.Restart).',1:6));
        assert(abs(max(A.PSValidated)-C.PSValidatedCase(i))<1e-9);
    end
    assert(abs(mean(C.SAValidated)-D.fourMethodTable.ValidatedAMI(4))<1e-9);
    assert(abs(mean(C.PSValidatedCase)-D.fourMethodTable.ValidatedAMI(3))<1e-9);
end

function check_synthetic_plot_export()
    T=synthetic_table();
    root=fso_result_utils.revision2_results_root();
    if exist(root,'dir')~=7,mkdir(root);end
    tmp=tempname(root);
    mkdir(tmp);
    c=onCleanup(@() remove_folder(tmp)); %#ok<NASGU>
    paths=plot_k5p2_fig7_four_methods(T,0.3,tmp,'CloseFigure',true);
    assert(exist(paths.pdf,'file')==2 && ...
        exist(paths.png,'file')==2 && exist(paths.fig,'file')==2);
end

function check_plot_rejects_bad_science()
    T=synthetic_table();
    root=fso_result_utils.revision2_results_root();
    tmp=tempname(root);
    bad=false;
    try
        plot_k5p2_fig7_four_methods(T,0.1,tmp,'CloseFigure',true);
    catch
        bad=true;
    end
    assert(bad,'Plot accepted wrong sigma_R^2.');
    T.Method(2)="SA";
    bad=false;
    try
        plot_k5p2_fig7_four_methods(T,0.3,tmp,'CloseFigure',true);
    catch
        bad=true;
    end
    assert(bad,'Plot accepted wrong method identities.');
    if exist(tmp,'dir')==7,remove_folder(tmp);end
end

function T=synthetic_table()
    % Deliberately artificial numbers; NEVER publication evidence.
    T=table(["Uniform PAM";"Zhou GS";"Pattern Search";"SA"], ...
        [0.85;1.06;1.15;1.25],[0;0;0.008;0.001], ...
        'VariableNames',{'Method','ValidatedAMI','StdAcrossReplicates'});
end
function remove_folder(path)
    if exist(path,'dir')==7,try,rmdir(path,'s');catch,end,end
end
