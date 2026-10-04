function F = i13d_benchmark_config()
%I13D_BENCHMARK_CONFIG  Frozen Reviewer-2 final optimizer benchmark.
%
% The benchmark compares the frozen I.11.1 SA results against MATLAB
% Pattern Search under a pre-declared, evaluation-budget-matched protocol.
%
% Frozen grid:
%   M          = [4 8 16 32]
%   SNR [dB]   = [20 30]
%   sigma_R^2  = [0 0.3]
%   replicates = [1 2]
%
% This yields 32 optimization cases. Each case uses the same number of
% starts as the frozen SA production: [6 6 8 8].
%
% Evaluation budget:
%   Pattern Search receives at most SA_maxIter+1 objective evaluations per
%   start (initial SA score + one proposal per SA iteration):
%       M=4,8,16 : 4001
%       M=32     : 8001
%
% Numerical equivalence threshold:
%   |I_SA-I_PS| <= 1e-3 bit/symbol.
%
% Scheduler:
%   dynamic parfeval queue, one restart per worker, no nested parallelism.
%   A new restart is dispatched immediately whenever one worker finishes.
%
% No global-optimum claim is implied.

    P=i11p1_production_config();

    F=struct();
    F.version='I13D-v1';
    F.description='Final SA-vs-Pattern-Search reviewer benchmark';
    F.MVec=[4 8 16 32];
    F.SNRVec=[20 30];
    F.TurbulenceVec=[0 0.3];
    F.Replicates=[1 2];
    F.P_avg=P.P_avg;
    F.minGap=P.minGap;
    F.pinZero=P.pinZero;
    F.tieToleranceBits=1e-3;
    F.ParallelWorkers=8;
    F.FreshPool=true;
    F.scheduler='parfeval-dynamic-v1';
    F.taskGranularity='one Pattern Search restart';
    F.noNestedParallelism=true;

    F.RestartsByM=zeros(size(F.MVec));
    F.MaxIterByM=zeros(size(F.MVec));
    F.MaxEvalsPerStart=zeros(size(F.MVec));
    for k=1:numel(F.MVec)
        [F.RestartsByM(k),~]=i11_restart_policy(F.MVec(k));
        [F.MaxIterByM(k),~]=i11p1_iteration_policy(F.MVec(k));
        F.MaxEvalsPerStart(k)=F.MaxIterByM(k)+1;
    end

    F.nPhysicalPoints=numel(F.MVec)*numel(F.SNRVec)*numel(F.TurbulenceVec);
    F.nOptimizationCases=F.nPhysicalPoints*numel(F.Replicates);
    F.nRestartTasks=0;
    F.totalMaxFunctionEvaluations=0;
    for k=1:numel(F.MVec)
        nCasesThisM=numel(F.SNRVec)*numel(F.TurbulenceVec)*numel(F.Replicates);
        F.nRestartTasks=F.nRestartTasks + nCasesThisM*F.RestartsByM(k);
        F.totalMaxFunctionEvaluations=F.totalMaxFunctionEvaluations + ...
            nCasesThisM*F.RestartsByM(k)*F.MaxEvalsPerStart(k);
    end

    % Pattern Search options. UseCompletePoll=false is deliberate: the
    % comparison is capped by an objective-evaluation budget, so a poll is
    % allowed to stop after finding an improving point rather than forcing
    % extra evaluations beyond the nominal cap.
    F.PatternSearch=struct();
    F.PatternSearch.Display='off';
    F.PatternSearch.UseCompletePoll=false;
    F.PatternSearch.MeshTolerance=1e-6;
    F.PatternSearch.StepTolerance=1e-8;
    F.PatternSearch.FunctionTolerance=1e-8;
    F.PatternSearch.UseParallel=false;

    % Longest-first dispatch policy. A task is prioritized first by whether
    % it is a turbulent refined-dt case, then by turbulent vs AWGN, then M,
    % SNR and sigma_R^2. Because dispatch is dynamic, this ordering actually
    % determines which single restart is sent next to the first free worker.
    F.priorityPolicy='refined-dt > turbulent > M > SNR > sigma_R^2';

    F.sourceProductionFreeze=P.freezeVersion;
    F.globalOptimumClaim=false;
end
