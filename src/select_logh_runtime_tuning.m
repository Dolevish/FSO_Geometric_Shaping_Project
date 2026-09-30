function [blockSize,symbolBatchSize,meta] = select_logh_runtime_tuning(M)
%SELECT_LOGH_RUNTIME_TUNING  Deterministic implementation-only tuning.
%
% Commit I.10.2.2 freezes the robust runtime policy measured by the
% I.10.2.1 production-grid tuning sweep.  The sweep covered representative
% low/medium/high SNR and turbulence points, two geometries, and all current
% production constellation orders M={4,8,16,32}.
%
% Selection rule: choose, for each M, the tested (block,batch) pair with the
% highest median speedup among candidates whose worst measured speedup was
% at least 0.98x versus the exact I.10.1 reference path.  This deliberately
% prefers robust settings over the more aggressive highest-median-only
% settings when the latter showed runtime regressions in some scenarios.
%
% These values affect only how identical floating-point likelihood work is
% grouped.  They do not change the y-grid, log-h quadrature, AMI definition,
% SA parameters, or feasible constellation set.

    validateattributes(M,{'numeric'},{'scalar','real','finite','integer','>=',2},mfilename,'M');
    M=double(M);

    if M <= 4
        blockSize = 2048;
        symbolBatchSize = 1;
    elseif M <= 8
        blockSize = 512;
        symbolBatchSize = 4;
    elseif M <= 16
        blockSize = 256;
        symbolBatchSize = 4;
    else
        % Current production maximum is M=32.  For larger future orders use
        % the measured M=32 robust setting conservatively until re-benchmarked.
        blockSize = 256;
        symbolBatchSize = 4;
    end

    symbolBatchSize=min(symbolBatchSize,M);

    meta=struct();
    meta.ruleVersion='I10.2.2-v1';
    meta.blockSize=blockSize;
    meta.symbolBatchSize=symbolBatchSize;
    meta.noRegressionThreshold=0.98;
    meta.source='I.10.2.1 measured production-grid robust tuning sweep';
end
