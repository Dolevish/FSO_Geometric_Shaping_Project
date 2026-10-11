# Commit K6 — Figure 8: 20 independent SA restart trajectories

## Purpose

Figure 8 (manuscript Section VI.E, "Repeatability of Simulated Annealing")
measures the variability of independently optimized SA trajectories under
ONE fixed channel condition. It is NOT a bar chart of 20 multi-start
replicate winners, and is not a new AMI-vs-SNR sweep.

Exactly 20 independently validated restart outputs, indexed 1 through 20:

- M=8
- SNR=20 dB
- sigma_R^2=0.1
- 4000 SA iterations PER restart
- frozen I.11.1 mean optical intensity 1, x(1)=0 and min gap 0.01
- equiprobable symbols, no instantaneous CSI, AWGN
- log-h fast evaluator with adaptive y grid, dt=0.01 for this condition
- independent adaptive integration validator for each restart

All SA settings remain identical to frozen I.11.1:

- T0=0.4, Tf=1e-3
- BaseStd0=0.15
- 50 iterations per temperature block
- same geometric cooling, dynamic proposal-step policy and full feasible
  projection, no new arbitrary xMax amplitude constraint
- history sampled every 50 iterations
- same sa_restart_task_checkpointed / sa_restart_task /
  simulated_annealing / build_fso_config and AMI_Validator.
- The original I.11.1 code, configs, MAT files and results are read-only.

The only deliberate design change relative to the per-replicate
I.11.1 execution is the number of restart trajectories (20 instead
of six) and a new predeclared deterministic random master seed.

### Independent restart RNG protocol

The K6 master seed is
uint32(fso_result_utils.case_seed(20261011,8,20,0.1)). The driver
explicitly rejects a seed that coincides with either source I.11.1
production replicate seed.

As in production sa_restart_task, restart 1 starts from uniform PAM.
Each restart r>=2 starts at PAM + 0.2 N(0,I), drawn from its individual
Threefry substream r. Its SA trajectory uses an independent Threefry
stream configured with the same master seed and substream r. The task's
outcome is independent of dynamic parallel scheduling order.

All 20 individual validated outcomes are included in the statistics.
There is NO picking a winner first, no discarding poor outcomes and no
averaging only multi-start winners.

## Statistics

For each restart r, I_r is the independently validated AMI of its
best-fast-AMI constellation. Compute:

    MeanAMI           = mean(I_1,...,I_20)
    SampleVarianceAMI = var(I_1,...,I_20,0)   [denominator 19]
    SampleStdAMI      = std(I_1,...,I_20,0)
    MinAMI            = min(I_r)
    MaxAMI            = max(I_r)
    RangeAMI          = MaxAMI - MinAMI
    MedianAMI         = median(I_r)

Additional diagnostics include worst/best restart, maximum absolute
fast-validator discrepancy, convergence K95/K99, runtime and acceptance
rates, all 8 intensity levels of each winner and actual adjacent spacing.

The K6 sample mean need not equal Table I/Fig. 7's SA mean, because
those summarize TWO multi-start validated winners instead of 20
individual random trajectory outcomes.

No global optimum is claimed.

## MATLAB runbook

Pull the updated ieee-revision branch in the repository root:

    git checkout ieee-revision
    git pull

In MATLAB with current folder set to the repository root:

    addpath(genpath(fullfile(pwd,'src')));

    test_k6_fig8_preflight();

    dry=sim_k6_fig8_repeatability('DryRun',true);

The preflight never runs optimization or evaluates AMI; the figure
export test uses clearly synthetic scores and discards its temporary
files. Every check should PASS before starting a real run.

Execute the production K6 trajectory experiment only after preflight:

    R6=sim_k6_fig8_repeatability();

Expected work:

- exactly 20 individually validated restart tasks;
- exactly 80,000 SA iterations;
- eight process workers by default, dynamic parfeval queue;
- every restart writes an atomic checkpoint, including full history.

On interruption, resume using the run directory printed in the
terminal; there is NO reoptimization of already validated checkpoints:

    R6=sim_k6_fig8_repeatability( ...
        'ResumeDirectory','<exact printed K6 run directory>');

If several I.11.1 production bundles exist, supply an explicit one:

    original=fullfile(fso_result_utils.default_results_root(), ...
      'i11p1_production','I11p1_production_20260930T224525960Z', ...
      'i11p1ProductionBundle.mat');

    test_k6_fig8_preflight('ProductionBundlePath',original);
    R6=sim_k6_fig8_repeatability('ProductionBundlePath',original);

The driver checks scientific configuration and original production
freeze against actual source metadata; mismatches fail before launch.
Resume requires the exact original signature and checkpoint identities.

## Result paths

All new files go under:

    results/ieee_revision2/commit_k6_fig8/K6_<UTC>/

- K6_repeatability_results.mat — complete metadata and each of 20
  detailed restart result structs, histories and constellations
- K6_restart_results.csv — one row per restart, seed, substream,
  fast/validated AMI, evaluation gap, elapsed times, convergence K95/K99,
  acceptance rate, actual minimum spacing and highest intensity
- K6_constellation_levels.csv — 20 rows with columns
  Restart, ValidatedAMI, x1, ..., x8
- K6_repeatability_summary.csv — one row for mean, sample variance,
  sample std, min, max, range, median and diagnostics
- Fig8_SA_repeatability.pdf — vector figure for IEEE manuscript
- Fig8_SA_repeatability.png — 400 DPI figure
- Fig8_SA_repeatability.fig — editable MATLAB figure
- K6_run_manifest.mat — exact frozen scientific signature
- checkpoints/K6_M8_SNR20_sig100_r01.mat through r20.mat

The bar chart has 20 bars in fixed restart-index order, a dashed
20-run sample-mean line and mean/std/variance/range annotations.
Its main y-axis starts at zero. When the observed spread is narrow,
an optional inset depicts that spread at a more legible scale without
altering the main zero-baseline chart.

### Manuscript wording

"Fig. 8. Validated AMI across 20 independent SA runs (M = 8,
SNR = 20 dB, sigma_R^2 = 0.1)."

Use measured K6 results to replace the [TBD] placeholders of
Section VI.E. Explicitly identify the sample variance (N-1)
and do not describe the spread as a global-optimality guarantee.

## Engineering/testing status

This commit changes only the K6 driver, figure generator, preflight and
this documentation. There is NO change to frozen I.11.1, K, K4 or K5
code or historical results. The GitHub code has not been executed
under MATLAB by the assistant. Run the real-source preflight first;
do not report synthetic plot-test values as experimental observations.
