# Commit K.3 — selected constellation four-panel figures

This post-processing step follows the visual structure of the project's legacy
constellation-comparison figure: intensity is shown on the x-axis and two stem/lollipop
series are overlaid in each panel.

Commit K has five d_min values. To obtain the requested four-panel layout, d_min=0 is
used as the unconstrained GS reference in every panel, while the four panels correspond
to:

    d_min = 0.005, 0.01, 0.02, 0.05.

A separate figure is generated for M=16 and M=32.

For every (M,d_min), including the d_min=0 reference, the plotted geometry is the
replicate having the highest validated AMI. No optimization or AMI evaluation is
rerun.

Each panel reports:

- imposed d_min,
- actual minimum adjacent spacing of the selected constellation,
- validated-AMI difference relative to the selected d_min=0 reference,
- reference and constrained validated AMIs.

All four panels of a given M share the same intensity-axis limits and ticks. This is
intentional: it keeps geometric movement directly comparable and avoids misleading
per-panel autoscaling.

Outputs are written by default to:

    <Commit-K run>/selected_constellation_figures/

Files:

    Fig_K_selected_constellations_M16.pdf
    Fig_K_selected_constellations_M16.png
    Fig_K_selected_constellations_M16.fig
    Fig_K_selected_constellations_M32.pdf
    Fig_K_selected_constellations_M32.png
    Fig_K_selected_constellations_M32.fig
    K_selected_constellation_plot_summary.csv

Run:

    Fk = plot_k_selected_constellations(Rk);

or:

    Fk = plot_k_selected_constellations('/path/to/K_min_gap_results.mat');

No global-optimum claim is made.
