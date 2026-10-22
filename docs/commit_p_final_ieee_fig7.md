# Commit P — Final IEEE Figure 7 (M8 SNR20 sigmaR2=0.3)

Commit P is a new, independent 12-start confirmation experiment at the same physical operating point as Commit M. It never overwrites M/N/O outputs.

## Frozen protocol
- M=8 equiprobable positive PAM, no-CSI IM/DD FSO log-normal fading and AWGN, SNR=20dB, sigma_R^2=0.3.
- Average optical intensity Pavg=1, x1=0, d_min=0.01.
- 12 distinct and feasible starts, with Uniform PAM first, 11 additional starts using independent fixed Threefry substreams and Gaussian perturbation std 0.2.
- For each restart r, SA and GPS receive the **same projected constellation**; the 12 starts are distinct from each other.
- The new, declared deterministic seed is `fso_result_utils.case_seed(P.baseSeed+104729,8,20,0.3)`, different from Commit M's seed. The offset 104729 is fixed in metadata and preflight tests. Cross-commit M/P restart indices are NOT paired.
- 12 SA tasks, each 4000 iterations and 4001 AMI evaluations counting the starting evaluation. Same cooling, proposal adjustments and Fast AMI as in previous experiments.
- 12 GPS tasks, **UseCompletePoll=true**, one fixed **MeshTolerance=1e-6**, StepTolerance=FunctionTolerance=1e-8, maximum 4001 objective evaluations. No additional mesh configurations and no GPS_CONTROL.
- Fast AMI log-h integration: dt=0.01, 8-sigma fading span, adaptive y-grid scale=1.1. Independently validated AMI for every final constellation.
- Uniform PAM and Zhou et al. GS are constructed anew at M=8, with Zhou ExtraBits=2 and Regularize=true; both are independently validated. They are not copied from Commit M.
- Four-method grayscale Figure 7 matches Commit O and M aesthetics, order Uniform PAM, Zhou GS, Pattern Search, SA. GPS and SA bars show **maximum validated AMI of all 12 restarts**. No error bars. Sample mean, variance, std, min and winner restart are saved separately across all 12 in CSV.
- Validate solver results, matching paired initial constellations, feasibility constraints and Fast/Validated gap. Save GPS exitflag, final mesh and objective function evaluations.
- Each actual task sends **START** and **DONE** with elapsed duration in minutes and seconds (e.g. `04m 32.16s`) via DataQueue. Baselines also log START and DONE.
- Load completed Commit M, N and O MAT files (automatic autodiscovery with optional `M8ResultPath`, `M16ResultPath`, `OResultPath` overrides). Save M/N/O/P comparisons and separate P-vs-M repeatability at equal M/SNR/sigmaR2. These comparisons do not determine which solution is an optimum; no global optimality claim.
- Use signed per-restart checkpoints and resume in place without rerunning finished tasks.

## MATLAB execution
```matlab
addpath(genpath(fullfile(pwd,'src')));
test_commit_p_final_ieee_preflight();
plan=run_commit_p_final_ieee_fig7('DryRun',true);
Pfinal=run_commit_p_final_ieee_fig7();
```
Optional override of previous result paths:
```matlab
Pfinal=run_commit_p_final_ieee_fig7('M8ResultPath','/path/CommitM_results.mat', ...
    'M16ResultPath','/path/CommitN_M16_results.mat', ...
    'OResultPath','/path/CommitO_results.mat');
```
For resuming: `'ResumeDirectory','/exact/results/ieee_revision2/commit_p_final_ieee_fig7/P_timestamp'`.
Optional `'UseParallel',false` for debugging; default is eight process workers.

## Files generated after MATLAB execution
```text
results/ieee_revision2/commit_p_final_ieee_fig7/P_<timestamp>/
  CommitP_results.mat
  CommitP_all_24_restarts.csv
  CommitP_optimizer_stats.csv
  Fig7_CommitP_four_methods.csv
  Fig7_CommitP_final.pdf
  Fig7_CommitP_final.png
  Fig7_CommitP_final.fig
  CommitP_comparison_M_N_O.csv
  CommitP_optimizer_comparison_M_N_O.csv
  CommitP_repeatability_vs_M.csv
  checkpoints/*
  native_gps/*
```

Commit P is a source-only commit, not an assertion of completed MATLAB numerical runs. Run on your MATLAB workstation with the Global Optimization and Parallel Computing Toolboxes.