# Commit K7 — standalone main and zoom figures for completed K6 Fig. 8

## Scope: FIGURE EXPORT ONLY

K6 has already completed 20 independently validated Simulated Annealing
restarts for manuscript Fig. 8, at M=8, SNR=20 dB, sigma_R^2=0.1.
K7 DOES NOT rerun SA, does not reevaluate AMI, and does not change the
original 20 scores, optimized constellations, scientific configuration,
I.11.1 production data, K6 checkpoints or historical figures.

The purpose is to resolve legibility issues in the original K6 figure:
the zoom inset obscured several bars and the statistics overlay covered
the upper portion of the data.

## Outputs: TWO separate figures

Main graph:
- Fig8_SA_repeatability_main.pdf
- Fig8_SA_repeatability_main.png
- Fig8_SA_repeatability_main.fig

Features:
- exactly 20 bars, in independent restart index order (1 to 20)
- main y axis starts at zero
- horizontal dashed line at the mean of the 20 validated AMIs
- the best restart has a slightly darker shade
- NO inset and no overlaid data annotations
- mean, sample std, sample variance and full range are printed in a
  dedicated band outside and above the plot axes

Zoom graph:
- Fig8_SA_repeatability_zoom.pdf
- Fig8_SA_repeatability_zoom.png
- Fig8_SA_repeatability_zoom.fig

Features:
- the same original 20 bars, same restart order and mean line
- standalone, full-size figure; NO inset
- zoomed y axis with a small padding outside the observed min/max
- explicitly marked "Expanded AMI axis (not zero-based)" to prevent
  a misleading interpretation of the magnitude of the differences

Both use MATLAB vector PDF export, 400-DPI PNG and editable FIG output.
No fake/reconstructed AMI values are used in publication exports.

### Statistical convention

The original K6 statistics are retained VERBATIM and verified against
the 20 actual measured AMI scores, including the sample (N-1) variance:

- mean = 1.706831690 bits/symbol
- sample variance = 1.28758334e-06 (bits/symbol)^2
- sample standard deviation = 0.00113471729 bits/symbol
- range = 0.00437932353 bits/symbol

These are from the user's completed K6 MATLAB stdout, and the exporter
checks the actual saved MAT source before creating any figure.
It does not use MATLAB population-variance settings (var(x,1)).
It never tries to extract "validatedAMI" as a property from the
K6 restartResults cell array; it reads the canonical K6 restartTable,
and crosschecks every saved restart result.

## Files added in K7

- src/export_k7_fig8_dual.m
- src/plot_k7_fig8_repeatability_dual.m
- src/test_k7_fig8_dual_export.m
- docs/commit_k7_fig8_dual_export.md

No existing project source file is edited, including
src/sim_k6_fig8_repeatability.m and src/plot_k6_fig8_repeatability.m.

## MATLAB commands

From the local repository root:

    git checkout ieee-revision
    git pull

Then in MATLAB:

    addpath(genpath(fullfile(pwd,'src')));

    test_k7_fig8_dual_export();

The preflight examines the completed K6 MAT and creates/discards
synthetic figure files only, to test both plots and vector export.
No SA, Pattern Search or AMI evaluation takes place.

For the completed run, use the still-existing MATLAB R6 workspace:

    K7 = export_k7_fig8_dual(R6);

Or load by explicit MAT path, even in a fresh MATLAB session:

    k6path = fullfile( ...
      fso_result_utils.revision2_results_root(), ...
      'commit_k6_fig8','K6_20261011T085135920Z', ...
      'K6_repeatability_results.mat');

    K7 = export_k7_fig8_dual(k6path);

Or auto-discover the most recently completed K6 run:

    K7 = export_k7_fig8_dual();

The exact source K6 directory used is printed. The saved MAT file is
only READ. No 80,000-iteration rerun is needed.

## Result location

K7 writes a NEW directory, separate from the original K6 result:

    results/ieee_revision2/
        commit_k7_fig8/
            K7_from_K6_20261011T085135920Z/
                Fig8_SA_repeatability_main.pdf
                Fig8_SA_repeatability_main.png
                Fig8_SA_repeatability_main.fig
                Fig8_SA_repeatability_zoom.pdf
                Fig8_SA_repeatability_zoom.png
                Fig8_SA_repeatability_zoom.fig
                K7_figure_manifest.mat

Using the same K6 source a second time replaces only the derived
K7 figure files and manifest; K6 original data/files remain untouched.

## Important limits

- K7 code was committed to GitHub, but needs MATLAB local preflight
  and visual inspection of the actual six outputs after the user runs it.
- A post-processing change does not alter the scientific interpretation:
  every bar remains one complete, independently validated 4000-iteration
  SA restart, not an average of multiple restarts.
- The larger zoom graph is intended to clarify visually small differences
  and must not be presented as a zero-based comparison.
- No global optimum claim.
