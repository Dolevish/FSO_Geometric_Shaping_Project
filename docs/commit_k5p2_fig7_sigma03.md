# Commit K5.2 — Frozen Fig. 7 for stronger log-normal turbulence

**Status:** Source code prepared for local MATLAB R2025b verification; no new MATLAB computation is claimed by this commit.

## Purpose

Reproduce the **same four-bar plot format as Commit K5** for
`M=8, SNR=20 dB, sigma_R^2=0.3`, rather than 0.1. All values come
from the completed I.13D *original* independently validated benchmark.
Default mode NEVER re-optimizes or evaluates AMI; it only audits and
aggregates frozen results.

Published-figure values are derived by re-selecting the validated
Pattern Search winner across all six restarts **for each replicate**, and
taking the two-replicate mean; the same procedure is used for the frozen
SA winners. Only PS and SA have between-replicate std error bars.
Zhou et al. GS and Uniform PAM are deterministic single-condition
baselines.

## Verified numerical reference (from previously exported I.13D CSV)

| Method | Mean validated AMI [bits/symbol] |
| --- | ---: |
| Uniform PAM | 0.9323308194 |
| Zhou GS | 1.1878905650 |
| Pattern Search | 1.2925719988 |
| SA | 1.3489855426 |

Expected SA minus PS ~0.05641354 bit/symbol at 0.3.
Commit K5 at 0.1 reported ~0.0101 bit/symbol (the exact ratio is read
from the saved K5 MAT, not hardcoded).

These are **descriptive**, previously observed results. sigma_R^2=0.3
was chosen after inspecting experiment outputs; avoid calling this
a pre-registered optimizer superiority test or a statistically established
advantage based on only 2 replicates.

## Scientific settings — UNCHANGED

Sources:
- I.11.1 production: saved validated SA, case seeds and PAM baselines.
- I.13B Zhou GS: saved 2025 literature GS evaluation.
- I.13D final benchmark: saved PS 224-task production, of which 12
  tasks at M8/20/0.3 are used.
- K5: saved two-replicate Fig. 7 at sigma_R^2=0.1 for comparison.

Assertions reject incompatible files or cases. All sources MUST point
to the same I.11.1 frozen production. Scientific policy remains:
- Average optical intensity (power) = 1, equiprobable symbols,
  fixed first level at 0, x_i>=0, adjacent level gap >= 0.01.
- 2 independent replicates, 6 restarts/replicate.
- 4000 frozen SA iterations/start, 4001 Pattern Search
  objective evaluations/start at most.
- Same raw Threefry start streams, exact I13D patternsearch worker,
  PS stopping tolerances and *direct linear feasible set*.
- Fast evaluator: log-h, candidate-adaptive y-grid, dt=0.01,
  YGridScale=1.1, span +/-8 sigma_t, block size 512.
- Validated restart selection (NOT fast-winner selection).
- Independent validator; no global-optimum claim.
- No changes to source implementations from I.11.1, I.13B,
  I.13D, K, K4 or K5.

The K5.2 primary figure is **always based on frozen I.13D** values.
An optional rerun of the exact PS tasks is an independent check and
is **never silently substituted** into Fig. 7.

## Files in this commit

- `src/make_k5p2_fig7.m`: default frozen replay, source auditing,
  re-aggregation, comparison against K5 and result export.
- `src/plot_k5p2_fig7_four_methods.m`: cloned K5 chart style, differing
  only in the 0.3 title and new output filename.
- `src/sim_k5p2_ps_recheck.m`: OPTIONAL same-policy 12-task
  PS recheck, dynamic 8-process `parfeval`, checkpoint/restart,
  START/DONE timing, independent validation and source comparison.
- `src/test_k5p2_fig7_preflight.m`: five structural/source/synthetic
  plot tests, none of which re-optimizes or evaluates AMI.

## Run in MATLAB (from repository root, branch ieee-revision)

```matlab
addpath(genpath(fullfile(pwd,'src')));
test_k5p2_fig7_preflight();
dry = make_k5p2_fig7('DryRun',true);
R52 = make_k5p2_fig7();
```

MAT files are on the user's computer, NOT tracked in GitHub.
Explicit paths are available if more than one candidate matches:
`'ProductionBundlePath'`, `'ZhouResultPath'`,
`'I13DResultPath'`, `'K5ResultPath'`.

For optional repeated PS validation **only if desired**:

```matlab
R52 = make_k5p2_fig7('RerunPatternSearch',true);
```

A stopped rerun can resume without redoing successful checkpoints:

```matlab
R52 = make_k5p2_fig7( ...
   'ResumeDirectory','<exact results/ieee_revision2/commit_k5p2_fig7/K5p2_...>', ...
   'RerunPatternSearch',true);
```

## Output files (results are NOT committed)

```text
results/ieee_revision2/commit_k5p2_fig7/K5p2_<timestamp>/
  Fig7_K5p2_four_methods.csv
  Fig7_K5p2_frozen_restarts.csv
  Fig7_K5p2_replicates.csv
  Fig7_K5p2_sigmaR2_0p3.pdf
  Fig7_K5p2_sigmaR2_0p3.png
  Fig7_K5p2_sigmaR2_0p3.fig
  K5_vs_K5p2_comparison.csv
  K5_vs_K5p2_SA_PS_gap.csv
  K5p2_results.mat
  ps_recheck/                           # only on optional rerun
    K5p2_PS_recheck_manifest.mat
    checkpoints/*.mat
    K5p2_rerun_patternsearch_restarts.csv
    K5p2_rerun_vs_I13D_replicates.csv
    K5p2_rerun_vs_I13D_summary.csv
    K5p2_rerun_results.mat
```

The chart is designed to retain the same four markers, spacing,
grayscale shading, axis label, standard-deviation error bars and
IEEE single-column style used for the verified K5 Figure 7.
Y-axis begins at zero.

**For submission:** If using 0.3 as Figure 7, update caption and
Section VI.D to 0.3 and avoid mixing the K5 (0.1) data in the
same numerical claim. This one selected condition is representative
and must not be presented as universal optimizer dominance.
