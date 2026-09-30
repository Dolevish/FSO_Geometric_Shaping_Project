function bundle = sim_i11_production_freeze(varargin)
%SIM_I11_PRODUCTION_FREEZE  Canonical final-production entry point.
%
% Scientific parameters are not exposed as name-value arguments here.  They
% come only from i11_production_config().  The caller may change operational
% controls (workers, output path, resume path) without changing the frozen
% scientific experiment.
%
% New run example:
%   b = sim_i11_production_freeze();
%
% Resume example (path printed at the start of the original run):
%   b = sim_i11_production_freeze('ResumeDirectory','/path/to/run');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'UseParallel',true,@logical_scalar);
    addParameter(p,'ParallelWorkers',8,@(v) isempty(v)||positive_integer(v));
    addParameter(p,'LongestFirst',true,@logical_scalar);
    addParameter(p,'WriteCSV',true,@logical_scalar);
    addParameter(p,'RunLabel','I11_production',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    addParameter(p,'ResumeDirectory','',@text_scalar);
    parse(p,varargin{:}); o=p.Results;

    F=i11_production_config();

    fprintf('\n============================================================\n');
    fprintf('I.11 PRODUCTION FREEZE (%s)\n',F.freezeVersion);
    fprintf('Grid: M=%s | SNR=%s dB | sigma_X^2=%s | replicates=%d\n', ...
        mat2str(F.MVec),mat2str(F.SNRVec),mat2str(F.TurbulenceVec),F.nReplicates);
    fprintf('SA: %d iterations | restarts by M=%s | T0=%.3g | Tf=%.3g | step0=%.3g\n', ...
        F.saMaxIter,mat2str(F.RestartsByM),F.T0,F.Tf,F.BaseStd0);
    fprintf('Expected optimization cases=%d | restart tasks=%d\n', ...
        F.nOptimizationCases,F.totalRestartTasks);
    fprintf('Scheduler: flat restart tasks | checkpoint/resume enabled | longest-first=%d\n',logical(o.LongestFirst));
    if o.UseParallel
        fprintf('Requested process workers: %s\n',mat2str(o.ParallelWorkers));
    end
    fprintf('No global-optimum claim is implied by this freeze.\n');
    fprintf('============================================================\n');

    % For the canonical local production run, honor an explicit worker count
    % even if an older pool with a different size is currently open.
    if o.UseParallel && ~isempty(o.ParallelWorkers)
        try
            pp=gcp('nocreate');
            if ~isempty(pp) && pp.NumWorkers~=double(o.ParallelWorkers)
                delete(pp);
            end
        catch
            % The underlying scheduler will handle pool creation/fallback.
        end
    end

    bundle=sim_revision_AMI_vs_SNR_checkpointed_flat( ...
        'MVec',F.MVec,'P_avg',F.P_avg,'SNRVec',F.SNRVec, ...
        'TurbulenceVec',F.TurbulenceVec,'minGap',F.minGap,'pinZero',F.pinZero, ...
        'FastFadingMethod',F.FastFadingMethod,'LoghDtPolicy',F.LoghDtPolicy, ...
        'loghDt',F.loghDt,'loghDtRefined',F.loghDtRefined, ...
        'loghSpanSigma',F.loghSpanSigma,'loghYBlockSize',F.loghYBlockSize, ...
        'YGridMode',F.YGridMode,'YGridScale',F.YGridScale, ...
        'saMaxIter',F.saMaxIter,'RestartPolicy',F.RestartPolicy, ...
        'nReplicates',F.nReplicates,'T0',F.T0,'Tf',F.Tf, ...
        'BaseStd0',F.BaseStd0,'ItersPerTemp',F.ItersPerTemp, ...
        'baseSeed',F.baseSeed,'historyEvery',F.historyEvery, ...
        'UseParallel',o.UseParallel,'ParallelWorkers',o.ParallelWorkers, ...
        'LongestFirst',o.LongestFirst,'Checkpoint',true, ...
        'ResumeDirectory',o.ResumeDirectory,'WriteCSV',o.WriteCSV, ...
        'RunLabel',o.RunLabel,'ResultsRoot',o.ResultsRoot);

    bundle.productionFreeze=F;
    bundle.productionFreeze.operational=struct('UseParallel',logical(o.UseParallel), ...
        'ParallelWorkers',o.ParallelWorkers,'LongestFirst',logical(o.LongestFirst), ...
        'RunLabel',char(string(o.RunLabel)));

    outPath=fullfile(bundle.run.outputDirectory,'i11ProductionBundle.mat');
    save(outPath,'bundle','-v7.3');
    fprintf('Saved canonical I.11 production bundle: %s\n',outPath);
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0&&mod(v,1)==0;end
