# Commit K5 — manuscript Fig. 7, frozen 0.1-turbulence four-method comparison

## Goal

Complete Figure 7, Section VI.D of the October 2026 manuscript draft,
for exactly M=8, SNR=20 dB, sigma_R^2=0.1, and compare validated AMI of:

1. Uniform PAM (frozen I.11.1 production baseline)
2. Zhou et al. literature GS (frozen I.13B baseline)
3. MATLAB constrained Pattern Search (the only new optimization)
4. Simulated Annealing (frozen I.11.1 production)

Existing final I.13D PS benchmarks sampled sigma_R^2={0,0.3}, so no final
budget-matched PS result was available at 0.1. The old short I.13C pilot
is NOT an acceptable replacement because its objective budget differs.

## Frozen scientific protocol — not changed

This commit adds new scripts ONLY; it does not edit any frozen I.11.1,
I.13B, I.13D, Commit K or Commit K4 source files or outputs.

- M=8; SNR=20 dB; sigma_R^2=0.1.
- P_avg=1, x_1=0 (pinZero), x_{i+1}-x_i >= d_min=0.01;
  equiprobable input levels, no additional artificial x_max constraint.
- Two original frozen production replicates, original production case seeds.
- Six restarts per replicate; raw initial points from the original
  i13_shared_start_points function (Threefry, same substreams as I.13D).
- I.13D patternsearch restart worker i13d_patternsearch_restart_task,
  not a new optimizer implementation.
- Pattern Search maximum 4001 fast-AMI function evaluations per restart
  (frozen M8 SA: 4000 iterations plus one initial score).
- Same fixed Pattern Search options and direct linear constraints:
  UseCompletePoll=false, UseParallel=false, MeshTolerance=1e-6,
  StepTolerance=1e-8, FunctionTolerance=1e-8,
  max final feasibility repair infinity norm=1e-5.
- Same log-h fast evaluator, candidate-adaptive y grid, YGridScale=1.1,
  8-sigma log-h span, y block size=512, frozen production dt=0.01.
- Only the canonical tiny post-termination projection repair allowed.
  Its post-hoc fast-score audit does not participate in the capped
  optimizer search budget.
- Each of 12 restart outputs is independently validated by the existing
  AMI_Validator, and each replicate winner is selected by validated AMI.
- For the four-method bar plot, SA and PS are **means of two independently
  validated replicate winners**, with sample std for those methods only.
  PAM and Zhou are deterministic at the chosen physical operating point.
- No claim of a global optimum, no new SA or Zhou optimizations.

The K5 driver asserts the frozen config values before running and checks
that the loaded I.13B Zhou and I.11.1 production files refer to the same
production run. Scientific settings and case identifiers are captured
in a K5 signature and checked on resume. Completed validated restarts
are checkpointed atomically by the existing I.13D task worker.

## Usage

In the updated repository on branch ieee-revision:

    git pull
    addpath(genpath(fullfile(pwd,'src')));

Run the preflight (no SA, Pattern Search or AMI evaluation):

    report = test_k5_fig7_preflight();

It checks the exact I.13D constraints/budgets, available frozen MATLAB
source data, two replicates with 6 restart keys each, unchanged Threefry
initialization and a synthetic-only plot export. It requires the actual
I.11.1 + I.13B MAT files on your computer.

Optional dry-run inspecting the resolved data and task plan:

    dry = fig7_m8_20db_sigma01('DryRun',true);

Then run ONLY the missing Pattern Search task grid:

    R5 = fig7_m8_20db_sigma01();

If the original source MAT paths are in an alternative location, provide
explicit ProductionBundlePath and ZhouResultPath for both the preflight
and the actual run; the source-production match remains compulsory.

On interruption, recover the exact saved K5 run directory displayed
in the log:

    R5 = fig7_m8_20db_sigma01( ...
        'ResumeDirectory','<exact K5 run directory>');

If a previous run successfully completed, the same command reuses all
12 validated checkpoints and regenerates only the tables/plot.

## Files generated

Default directory:

    results/ieee_revision2/commit_k5_fig7/K5_<UTC>/

- Fig7_patternsearch_restart_results.csv
- Fig7_four_method_values.csv
- Fig7_four_method_results.mat
- Fig7_four_methods_M8_20dB_sigmaR2_0p1.pdf
- Fig7_four_methods_M8_20dB_sigmaR2_0p1.png
- Fig7_four_methods_M8_20dB_sigmaR2_0p1.fig
- K5_scientific_manifest.mat
- checkpoints/*.mat (12)

The graph is for one physical case, so bars are more appropriate than
SNR- or turbulence-curves. Results from sigma_R^2=0.3 are intentionally
NOT mixed into this figure.

## Status and manuscript usage

The code has been committed, but no MATLAB production experiment was
executed or certified during this code change. Only after the local
MATLAB preflight and all 12 validated PS restarts finish should values
from Fig7_four_method_values.csv be inserted in Section VI.D and its
[ TBD ] entries updated. The manuscript's caption remains
"Validated AMI for uniform PAM, literature GS, Pattern Search and SA
(M=8, 20 dB, sigma_R^2=0.1)." Verify final figure legibility at
IEEE single-column scale.
