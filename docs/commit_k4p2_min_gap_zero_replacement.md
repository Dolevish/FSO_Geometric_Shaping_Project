# Commit K.4.2 — replace zero-minimum-spacing with measured 0.001

## Frozen scientific setup

Add the two **new** physical points `M=[16,32]`,
`SNR=20 dB`, `sigma_R^2=0.2`, `d_min=0.001`.
The same `d_min=0.005` measurements **already exist** in Commit K,
so they are deliberately NOT rerun. Existing K / K4 data are read-only.

All scientific/numerical controls are inherited without alteration from
`k_min_gap_config()` and the original K runner:

- \(M=16,32\), \(SNR=20\) dB, \(\sigma_R^2=0.2\).
- Fixed **mean optical intensity** \(P_{\mathrm{avg}}=1\) (not squared
  electrical energy), \(x_1=0\), ordered nonnegative constellation levels.
- 2 independent replicates using exactly the original paired
  `fso_result_utils.case_seed` seeds and Threefry restart streams.
- 8 **independent restarts per replicate**, **8000 SA iterations per start**.
- Temperature profile `T0`, `Tf`, `BaseStd0`,
  `ItersPerTemp`, same cooling and proposal mechanisms.
- Accurate original log-h fast-AMI evaluator, candidate-adaptive y grid,
  frozen selected `dt=0.01`, independent original AMI validator.
- Every restart independently validated. Replicate winner determined
  by maximum **validated** AMI, not fast-evaluator AMI.
- The original 8-worker, dynamic `parfeval` queue, START/DONE logging,
  task/SA/validation durations, atomic checkpoints, resume support.
- No global-optimum claim.

The only change to an existing source file is **two narrow gateway
extensions** in `sim_k_min_gap_sensitivity.m`: allow the new
`K4p2-v1` FrozenConfig version and use a unique `K4p2_extension`
output stem. Original K/K4 branches retain precisely their behavior.

## New point count

| M | new d_min | replicates | restarts/rep | task count |
|---|---|---:|---:|---:|
| 16 | 0.001 | 2 | 8 | 16 |
| 32 | 0.001 | 2 | 8 | 16 |

Total 2 physical points, 4 replicate cases, **32 restarts**,
**256,000 SA iterations** and **256,032 nominal fast evaluations**.

## New merged publication dataset

Take the **existing 16-point K4 combined results**, remove the
zero-gap rows **only from the new derived dataset**, and add the
2 validated d_min=0.001 rows. The 0.005 rows are reused unchanged.
Final point counts: 9 at M16; 7 at M32. Total 16 physical
points, 32 replicate cases and 256 restart rows.

M=16 merged d_min:
`[.001,.005,.01,.02,.03,.04,.05,.07,.10]`

M=32 merged d_min:
`[.001,.005,.01,.02,.03,.04,.05]`

The old metric `GainRetentionVsGap0` is inappropriate once d_min=0
is excluded; the new dataset replaces it with an explicitly named
`GainRetentionVsGap0001`. Existing mean validated AMI and replicate
std of the historical nonzero points are copied exactly without
re-evaluation. The original source MAT/CSV files are not deleted or
modified.

## Plot style

`plot_k4p2_min_gap` follows the original K4 combined plot:
- one axis, both M curves, same blue circles / orange squares,
  Times New Roman, labeled axes, legend, grid, line weights,
  feasibility shading beyond M32's maximum possible gap;
- no error bars (sample standard deviations remain in the CSV);
- exactly the new values including 0.001; no zero point or fake
  M32 points at 0.07/0.1.
- exported as PDF, PNG 400dpi, and editable MATLAB FIG.

**Exact custom-source FIG option:** since MATLAB FIG contains
application-specific graphics objects, the renderer additionally
supports using your own existing
`Fig_K4_combined_min_gap_no_errorbars(1).fig` as a visual template.
Supply `'TemplateFigPath','/path/to/original.fig'`; it opens the
FIG, substitutes the two curves' XData/YData only, updates axis
limits/ticks, and saves **new** output files. The original figure is
never overwritten. Without this option the code reproduces the K4
plot styling directly.

## Setup and preflight

```bash
git checkout ieee-revision
git pull
```

In MATLAB from project root:

```matlab
addpath(genpath(fullfile(pwd,'src')));
test_k_min_gap_sensitivity();
test_k4_min_gap_extension();
test_k4p2_min_gap();
D42 = sim_k4p2_min_gap('DryRun',true);
```

Preflight uses a **synthetic added-case fixture copied from the
historical 0.005 source row** to audit the merged result and PNG/PDF/FIG
export; these are NOT new optimization results and are deleted after
the test. This preflight requires the original frozen K4 merged MAT.
If discovery is ambiguous:

```matlab
test_k4p2_min_gap('K4CombinedPath',originalK4Mat);
```

After passing preflight run ONLY new 0.001 cases:

```matlab
R42 = sim_k4p2_min_gap();
```

To resume without recomputing completed validated restarts:

```matlab
R42 = sim_k4p2_min_gap( ...
    'ResumeDirectory','<actual K4p2 K_... run directory>');
```

## Merge and export

Use the previously completed K4 combined MAT:

```matlab
root=fso_result_utils.revision2_results_root();
originalK4Mat=fullfile(root,'commit_k_min_gap_extension', ...
    'K_20261009T081312445Z','combined','K4_combined_results.mat');
C42=merge_k4p2_min_gap(originalK4Mat,R42);
```

Optionally apply your exact previous no-errorbars MATLAB FIG template:

```matlab
C42=merge_k4p2_min_gap(originalK4Mat,R42, ...
    'TemplateFigPath','/absolute/path/Fig_K4_combined_min_gap_no_errorbars(1).fig');
```

The second command writes new outputs in the K4p2 run directory;
it never changes the template FIG or old MAT. The source K4 plot
is not loaded unless explicitly requested.

## Output paths

```text
results/ieee_revision2/commit_k4p2_min_gap/K_<timestamp>/
  K4p2_extension_results.mat
  K4p2_extension_restart_results.csv
  K4p2_extension_replicate_summary.csv
  K4p2_extension_physical_summary.csv
  checkpoints/
  combined/
    K4p2_combined_results.mat
    K4p2_combined_restart_results.csv
    K4p2_combined_replicate_summary.csv
    K4p2_combined_physical_summary.csv
    Fig_K4p2_combined_min_gap_no_errorbars.pdf
    Fig_K4p2_combined_min_gap_no_errorbars.png
    Fig_K4p2_combined_min_gap_no_errorbars.fig
```

No production optimization has been executed as part of the Git
commit; run the MATLAB preflight and inspect the output before
reusing it in the reviewer response or revised paper.
