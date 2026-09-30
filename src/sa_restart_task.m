function result = sa_restart_task(cfg, x_baseline, masterSeed, restartIndex, taskLabel)
%SA_RESTART_TASK  Run one deterministic SA restart plus independent validation.
%
% Commit I.10.5A extracts the atomic restart unit from sa_multistart so
% restart tasks from different physical cases can share one global parallel
% queue.  The random-number construction deliberately matches sa_multistart:
%
%   - the case master seed is converted to uint32;
%   - restart 1 starts from the supplied baseline;
%   - restart k>=2 uses Threefry substream k to perturb the baseline;
%   - the SA trajectory itself uses a fresh Threefry stream with the SAME
%     master seed and substream k.
%
% Therefore the result of a restart is independent of execution order and
% can be reproduced whether it runs serially, inside sa_multistart, or on a
% process worker in a flattened restart-task scheduler.
%
% The returned structure is field-compatible with one element of the
% sa_multistart results array, with additional task-timing fields.

    if nargin < 5 || isempty(taskLabel)
        taskLabel = '[SA]';
    end
    if isstring(taskLabel), taskLabel = char(taskLabel); end

    validateattributes(restartIndex,{'numeric'}, ...
        {'scalar','integer','positive','finite'},mfilename,'restartIndex');
    if ~isfield(cfg,'AMI_Validator') || isempty(cfg.AMI_Validator)
        error('sa_restart_task:MissingValidator', ...
            'cfg.AMI_Validator is required for validated restart execution.');
    end

    masterSeed = uint32(masterSeed);
    restartIndex = double(restartIndex);
    x_baseline = x_baseline(:);

    if restartIndex == 1
        x0_raw = x_baseline;
    else
        sX0 = RandStream('Threefry','Seed',double(masterSeed));
        sX0.Substream = restartIndex;
        x0_raw = x_baseline + 0.2*randn(sX0,numel(x_baseline),1);
    end

    % Match sa_multistart exactly: a separate Threefry stream drives SA.
    s = RandStream('Threefry','Seed',double(masterSeed));
    s.Substream = restartIndex;
    RandStream.setGlobalStream(s);

    label = sprintf('%s-%02d',taskLabel,restartIndex);
    taskStartPosix = posixtime(datetime('now','TimeZone','UTC'));
    tTask = tic;
    outk = simulated_annealing(cfg,x0_raw,restartIndex,label);

    AMI_functions.assert_constellation_feasible(outk.x_best,cfg,1e-10);
    tVal = tic;
    miVal = cfg.AMI_Validator(outk.x_best);
    validationRuntime = toc(tVal);
    if ~(isscalar(miVal) && isfinite(miVal))
        error('sa_restart_task:InvalidValidatedAMI', ...
            'Validation returned an invalid AMI for restart %d.',restartIndex);
    end

    result = struct();
    result.startIndex = restartIndex;
    result.masterSeed = masterSeed;
    result.substream = restartIndex;
    result.x0_raw = x0_raw(:);
    result.x0_projected = outk.x0_projected(:);
    result.initialMI = outk.initialMI;
    result.bestMI = outk.bestMI;
    result.bestMIFast = outk.bestMI;
    result.bestMIValidated = miVal;
    result.validationGap = outk.bestMI-miVal;
    result.validationRuntime = validationRuntime;
    result.fastRank = NaN;
    result.validatedRank = NaN;
    result.x_best = outk.x_best(:);
    result.finalMI = outk.finalMI;
    result.finalX = outk.finalX(:);
    result.runtime = outk.runtime;
    result.nIterations = outk.nIterations;
    result.nAccepted = outk.nAccepted;
    result.nImprovingAccepted = outk.nImprovingAccepted;
    result.nWorseAccepted = outk.nWorseAccepted;
    result.nRejected = outk.nRejected;
    result.acceptanceRate = outk.acceptanceRate;
    result.history = outk.history;
    result.taskWallSeconds = toc(tTask);
    result.taskStartPosix = taskStartPosix;
    result.taskEndPosix = posixtime(datetime('now','TimeZone','UTC'));

    fprintf('%s validate | fast=%.6f | val=%.6f | delta=%+.3e | %.2fs\n', ...
        label,result.bestMIFast,result.bestMIValidated,result.validationGap,validationRuntime);
end
