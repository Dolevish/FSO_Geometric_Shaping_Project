# Commit K8 — Best-validated AMI curves in original IEEE visual style

## Scope and scientific definition

Change **only the plotted performance-point selection**, from the historical
mean of the two I.11.1 production replicate winners to the **highest observed
validated AMI among ALL saved SA restarts** for each physical channel point
`(M, SNR_dB, sigma_R^2)`.

The canonical frozen I.11.1 run covers:
- M = 4, 8, 16, 32; SNR = 5:5:30 dB; sigma_R^2 = 0, 0.1, 0.2, 0.3;
- two independent production replicates at every physical point;
- six SA restarts per replicate for M=4,8 and eight for M=16,32;
- 4000 SA iterations per restart for M=4,8,16; 8000 for M=32;
- winner **within each replicate selected by independent validated AMI**.

Consequently, for each point the maximum over ALL validated restarts equals
the maximum of the two **already-validated** replicate winners. The new
curve is an **observed best-run** curve, not an average or a global optimum.

```
I_best(M,SNR,sigma) = max( I_validated(rep1), I_validated(rep2) )
Gain_best = I_best - I_PAM_validated
```

Tie policy: deterministic smallest replicate index (replicate 1).

The Uniform PAM reference curve stays identical. The noise/fading model,
normalizations, M/SNR/turbulence grid, SA, validation code, seed schedule,
and all frozen production MAT files remain unchanged.

## Exact visual formatting

The existing function `make_ieee_revision_figures_legacy_style.m` now
supports `PerformanceAggregation='best-validated'`. Both original paper
performance figure renderers and styles are reused:
- **AMI vs SNR** (M=8, 4 turbulence curves): colors, line and marker
  specifications, legend labels/position, fonts, ticks, limits unchanged.
- **AMI vs SNR allM** (4 turbulence subplots and all M): same markers,
  colors, layout, axes, legend, line weights and automatic upper bound policy.
- Vector PDF, 300 DPI PNG, and editable MATLAB FIG retain the **same
  figure filenames** but live under a **different folder**.

The separate constellation-comparison figure (Fig. 4) is not regenerated
by the best-curve mode: its original representative-replicate-1 content
is unrelated to selecting a winner at every curve point. Mean mode still
exports all three original figures exactly as before.

## How to run from repository root, on ieee-revision

```matlab
git checkout ieee-revision % run in shell, NOT inside MATLAB
git pull              % run in shell
```

In MATLAB:

```matlab
addpath(genpath(fullfile(pwd,'src')));
test_best_validated_curve_selection(); % synthetic-only, no SA or AMI
out = make_ieee_revision_best_validated_legacy_curves('CloseFigures',true);
disp(out.amiVsSNR.pdf);
disp(out.amiVsSNRAllM.pdf);
disp(out.bestValidatedTableCSV);
```

If auto-discovery is ambiguous, supply the **actual frozen production MAT**
by full path using `'ProductionBundlePath','.../i11p1ProductionBundle.mat'`.

An alternative equivalent command is:

```matlab
out = make_ieee_revision_figures_legacy_style( ...
    'PerformanceAggregation','best-validated','CloseFigures',true);
```

## Outputs

Default new path:
`<production>/i12_analysis/ieee_figures_legacy_style_best_validated/`

- `AMI vs SNR.pdf`, `AMI vs SNR.png`, `AMI vs SNR.fig`
- `AMI vs SNR allM.pdf`, `AMI vs SNR allM.png`, `AMI vs SNR allM.fig`
- `best_validated_curve_points.csv` (96 rows, PAM, each replicate validated
  AMI, absolute winner, winner replicate, gain, and restarts considered)
- `i12c_legacy_figure_manifest.mat` recording exact aggregation and source.

Original mean-performance figures:
`<production>/i12_analysis/ieee_figures_legacy_style/` remain untouched.
The original mean setting also remains callable by the old no-arg API.

## Validation and publication caveat

- The selector checks for all 192 frozen case rows, exactly both replicate
  identities per physical point, valid finite AMI, consistent PAM reference
  and agreement with the previously analyzed I.12A mean/replicate columns.
- The source case table already aggregates all **validated** restart winners.
  No optimization or AMI function is ever called by this change.
- `test_best_validated_curve_selection` checks numeric selection, tie,
  permutation invariance and failures on incomplete/inconsistent input.
- Actual graphs can only be exported on the machine holding the frozen
  `i11p1ProductionBundle.mat`, because generated MATLAB results are
  ignored by Git and not checked into the repository.
- When quoting a graph, state **best of 12 or 16 SA restarts** for each
  physical point, not mean of replicates or global optimum.
