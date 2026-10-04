# I.13E-F — Reviewer 1/2 publication outputs

This stage is post-processing only. It does not run SA, Pattern Search, the fast AMI
objective, or the independent validator.

## Purpose

I.13E-F converts the frozen I.13B/I.13D results into publication-ready material for
Reviewer 1 and Reviewer 2.

Reviewer 1 is addressed by separating:

1. generic geometric-shaping gain relative to uniform PAM, using the deterministic
   Zhou-2025 optical-intensity GS literature baseline; and
2. additional channel-adaptation gain obtained by optimizing directly for the
   marginalized no-instantaneous-CSI log-normal FSO AMI.

Reviewer 2 is addressed by comparing the frozen SA solution with constrained Pattern
Search under the I.13D predeclared protocol.

No global-optimum claim is made.

## Figures

### Fig_I13EF_optimizer_comparison_ieee

Two-panel white-background IEEE-style figure.

- Panel (a): validated Pattern Search AMI versus validated SA AMI, with y=x.
- Panel (b): I_SA-I_PS grouped by M and by the four benchmark channel conditions.
  Error bars show the sample standard deviation of the two replicate deltas.
  The +/-1e-3 bit/symbol predeclared equivalence region is shaded.

### Fig_I13EF_four_method_ieee

Four panels, one per M, comparing:

- uniform PAM;
- Zhou literature GS;
- channel-adaptive Pattern Search GS;
- channel-adaptive SA GS.

The four operating conditions are shown as categorical groups. No lines connect
20/0, 20/0.3, 30/0 and 30/0.3, avoiding the false visual implication that this
categorical x-axis is continuous. Error bars are shown for the two stochastic
optimizers using their two-replicate sample standard deviations.

## Tables

- Table_I13EF_all_pointwise.csv
  Complete 16-point comparison with PAM, Zhou, Pattern Search, SA, all gain
  decompositions and optimizer outcome.

- Table_I13EF_optimizer_summary_by_M.csv
  Reviewer-2 summary by constellation order, including outcome counts, delta
  statistics and restart/replicate variability.

- Table_I13EF_turbulent_method_comparison.csv
  The eight sigma_R^2=0.3 operating points. This is the compact candidate table for
  a response letter or manuscript supplement.

- Table_I13EF_turbulent_gain_summary_by_M.csv
  Mean turbulent gains by M:
    Zhou-PAM,
    PatternSearch-Zhou,
    SA-Zhou,
    SA-PatternSearch.

- I13EF_reviewer12_summary.txt
  Automatically generated numerical narrative/checklist. It explicitly distinguishes
  the identical maximum evaluation budget from the actual number of Pattern Search
  search evaluations, since Pattern Search may stop early.

## Recommended call

    O = make_i13ef_reviewer12_publication_outputs(Rfinal);

or, after reopening MATLAB,

    O = make_i13ef_reviewer12_publication_outputs( ...
        '/path/to/i13d_optimizer_benchmark.mat');

By default all outputs are placed in:

    <I13D run directory>/publication_ready
