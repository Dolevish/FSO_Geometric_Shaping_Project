function R42 = sim_k4p2_min_gap(varargin)
%SIM_K4P2_MIN_GAP  Run ONLY M16/M32 at d_min=0.001.
%
% Uses the original frozen K simulation runner (same SA, evaluator,
% validator, Threefry restarts, checkpoint protocol, dynamic 8-worker pool).
% Existing K/K4 results and d_min=0.005 checkpoints are not touched.
%
%   dry = sim_k4p2_min_gap('DryRun',true);
%   R42 = sim_k4p2_min_gap();
%   R42 = sim_k4p2_min_gap('ResumeDirectory','/path/to/K_...');

    K=k4p2_min_gap_config();
    parent=fullfile(fso_result_utils.revision2_results_root(), ...
        'commit_k4p2_min_gap');
    R42=sim_k_min_gap_sensitivity( ...
        'FrozenConfig',K, ...
        'OutputDirectory',parent, ...
        'MakeFigure',false, ...
        varargin{:});
end
