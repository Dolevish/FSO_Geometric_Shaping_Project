function [blockSize,symbolBatchSize,meta] = select_logh_runtime_tuning(M)
%SELECT_LOGH_RUNTIME_TUNING  Deterministic implementation-only tuning.
%
% Commit I.10.2 uses the I.10.1 runtime benchmark to choose a conservative
% y-block size by constellation order, plus a small symbol batch for the
% likelihood kernel.  These values change only how the same floating-point
% work is grouped; they do not change the y-grid, log-h quadrature, SA, or
% feasible constellation set.
%
% Current production orders are M={4,8,16,32}.  The block-size rule is based
% on the measured I.10.1 best blocks for M=4,16,32; M=8 inherits the small-M
% setting.  Symbol batching is deliberately conservative to control memory.

    validateattributes(M,{'numeric'},{'scalar','real','finite','integer','>=',2},mfilename,'M');
    M=double(M);

    if M <= 8
        blockSize = 1024;
    elseif M <= 16
        blockSize = 2048;
    else
        blockSize = 256;
    end

    if M <= 8
        symbolBatchSize = 2;
    elseif M <= 16
        symbolBatchSize = 1;
    else
        symbolBatchSize = 2;
    end
    symbolBatchSize=min(symbolBatchSize,M);

    meta=struct();
    meta.ruleVersion='I10.2-v1';
    meta.blockSize=blockSize;
    meta.symbolBatchSize=symbolBatchSize;
    meta.source='I10.1 measured block-size benchmark plus conservative memory-aware batching';
end
