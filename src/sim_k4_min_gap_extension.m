function R4 = sim_k4_min_gap_extension(varargin)
%SIM_K4_MIN_GAP_EXTENSION Run only the six feasible new spacing points.
%
% Scientific settings are inherited from frozen Commit K. Uses the same
% checkpointed dynamic scheduler, SA starts, objective, and validator.
%
%   R4=sim_k4_min_gap_extension('DryRun',true);
%   R4=sim_k4_min_gap_extension();
%   R4=sim_k4_min_gap_extension('ResumeDirectory','/path/to/K_...');
%
% Outputs under results/ieee_revision2/commit_k_min_gap_extension/K_<UTC>.
% No rerun of the 10 original Commit-K physical points is performed.

    K=k4_min_gap_extension_config();
    outputParent=fullfile(fso_result_utils.revision2_results_root(), ...
        'commit_k_min_gap_extension');

    % Explicit caller options follow these defaults (inputParser keeps the
    % standard ResumeDirectory, ParallelWorkers, DryRun etc.).
    R4=sim_k_min_gap_sensitivity( ...
        'FrozenConfig',K, ...
        'OutputDirectory',outputParent, ...
        'MakeFigure',false, ...
        varargin{:});
end
