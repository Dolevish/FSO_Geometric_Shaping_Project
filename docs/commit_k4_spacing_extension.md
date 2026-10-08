# Commit K.4 — six-point spacing extension and merged figure

## Decisions and scientific scope

The old completed Commit K data are immutable. The extension runs only six
additional feasible (M,d_min) physical points:

| M | new d_min |
|---|---|
| 16 | 0.03, 0.04, 0.07, 0.10 |
| 32 | 0.03, 0.04 |

M=32 cannot use d_min=0.07 or 0.10 under the unchanged mean optical-power
constraint P_avg=1 with the first level pinned to zero. Necessary feasibility:

    d_min <= 2 * P_avg / (M-1)

All other scientific inputs are inherited unchanged from frozen Commit K:
SNR=20dB, sigma_R^2=0.2, P_avg=1, two replicates, eight restarts per
replicate, 8000 iterations/restart, SA cooling and proposal policy, same
candidate-adaptive fast likelihood evaluator and independent validator.

The task budget is exactly 6 physical points, 12 replicate cases,
96 validated restarts, 768,000 SA iterations, 768,096 nominal fast
objective evaluations. Eight-process dynamic queue, START/DONE timing logs,
atomic checkpoints, deterministic paired seeds and exact signature matching
are the same as Commit K.

The old K-v1 config and signature remain unchanged, and no existing run is
resumed into a scientifically different set of tasks.

## 1. Pull and preflight

From the repository root:

    git pull

In MATLAB:

    addpath(genpath(fullfile(pwd,'src')));
    test_k_min_gap_sensitivity();   % original frozen K regression
    test_k4_min_gap_extension();    % 6 K4 checks, includes synthetic merge

No SA is performed by these structural tests.

Inspect the plan without running:

    dry=sim_k4_min_gap_extension('DryRun',true);

## 2. Execute the six new physical points

    R4=sim_k4_min_gap_extension();

This is the only computationally expensive command. It dispatches exactly
96 new validated restart tasks with the same dynamic scheduler. Expected
duration depends on actual hardware and workload. No original task is rerun.

On interruption use the actual run directory printed in the MATLAB log:

    R4=sim_k4_min_gap_extension( ...
        'ResumeDirectory','<printed K4 run directory>');

Extension results are saved at:

    results/ieee_revision2/commit_k_min_gap_extension/K_<UTC>/

Important files:

    K4_extension_results.mat
    K4_extension_restart_results.csv
    K4_extension_replicate_summary.csv
    K4_extension_physical_summary.csv
    K_scheduler_trace.csv
    K_run_manifest.mat
    checkpoints/

## 3. Merge with the completed original Commit K result

The original run is:

    results/ieee_revision2/commit_k_min_gap/K_20261006T084906791Z/

Use its MAT-file path, without rerunning that experiment:

    originalK=fullfile( ...
      fso_result_utils.revision2_results_root(), ...
      'commit_k_min_gap','K_20261006T084906791Z','K_min_gap_results.mat');

    C4=merge_k4_min_gap_results(originalK,R4);

The merger verifies scientific settings, baseline PAM values, paired seeds,
unique case/restart identities, all validation winners, fixed-power and
spacing constraints, and that the original table rows are unchanged.

New gain-retention statistics are computed relative to the original
d_min=0 point for each M.

After merging the complete dataset comprises:

- M=16: 9 points
- M=32: 7 points
- Total: 16 physical points, 32 replicate cases, 256 restart tasks

The merged MAT/CSV files are saved in:

    <K4 run>/combined/

with names:

    K4_combined_results.mat
    K4_combined_restart_results.csv
    K4_combined_replicate_summary.csv
    K4_combined_physical_summary.csv

## 4. Figure

The merger automatically calls plot_k4_combined_min_gap and creates:

    Fig_K4_combined_min_gap.pdf
    Fig_K4_combined_min_gap.png
    Fig_K4_combined_min_gap.fig

The figure shows mean validated AMI across two replicate winners with
sample-standard-deviation error bars. M=32 stops at d_min=0.05. A light
region and vertical boundary distinguish M=32's infeasible domain beyond
d_min=2/31. The figure does not extrapolate unknown results.

Both original and extension MAT files remain unmodified. Neither SA nor AMI
is rerun during merge/plotting. No global-optimum claim is made.
