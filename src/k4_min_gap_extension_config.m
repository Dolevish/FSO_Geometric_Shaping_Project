function K = k4_min_gap_extension_config()
%K4_MIN_GAP_EXTENSION_CONFIG  Frozen six feasible extension points.
%
% M=16: d_min=[0.03 0.04 0.07 0.10]
% M=32: d_min=[0.03 0.04] (larger gaps violate P_avg=1).
%
% ALL scientific parameters are inherited unchanged from Commit K.
% The original K-v1 configuration and completed results stay untouched.

    K=k_min_gap_config();
    K.version='K4-v1';
    K.description='Commit K extension: six additional feasible minimum-spacing points';
    K.MinGapVec=[0.03 0.04 0.07 0.10];
    K.MinGapByM={[0.03 0.04 0.07 0.10],[0.03 0.04]};

    if ~isequal(K.MVec,[16 32])
        error('k4_min_gap_extension_config:M','Unexpected M grid.');
    end
    nPhysical=0;
    for i=1:numel(K.MVec)
        m=K.MVec(i);
        gaps=K.MinGapByM{i};
        if any(gaps<0)||any(gaps>2*K.P_avg/(m-1)+1e-12)
            error('k4_min_gap_extension_config:Infeasible', ...
                'Infeasible d_min for M=%d, P_avg=%g.',m,K.P_avg);
        end
        nPhysical=nPhysical+numel(gaps);
    end

    K.nPhysicalSensitivityPoints=nPhysical;
    K.nReplicateCases=nPhysical*numel(K.Replicates);
    K.nRestartTasks=K.nReplicateCases*K.NStarts;
    K.totalSAIterations=K.nRestartTasks*K.MaxIter;
    K.nominalFastObjectiveEvaluations=K.nRestartTasks*(K.MaxIter+1);
    K.globalOptimumClaim=false;

    assert(K.nPhysicalSensitivityPoints==6 && K.nReplicateCases==12 && ...
        K.nRestartTasks==96 && K.totalSAIterations==768000);
end
