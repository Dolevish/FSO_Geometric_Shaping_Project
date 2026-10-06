# Commit K.2 — selected-constellation reporting

This post-processing step uses the already completed Commit-K result. It does not rerun
SA, the fast AMI evaluator, or the independent validator.

For every physical point (M,d_min), the replicate with the highest validated AMI is
selected. Ties are broken deterministically by choosing the lower replicate number.

Run:

    Qk = report_k_selected_constellations(Rk);

or after restarting MATLAB:

    Qk = report_k_selected_constellations( ...
        '/path/to/K_min_gap_results.mat');

The function prints three tables:

1. Final constellation levels for each selected (M,d_min) case.
2. Compact spacing table with M, imposed d_min, actual d_min, and validated AMI.
3. Detailed table with M, SNR, sigma_R^2, selected replicate/restart, validated and
   fast AMI, actual minimum spacing, slack, binding flag, number of adjacent pairs
   on the imposed spacing, maximum level, mean/median adjacent gap, complete levels,
   and complete adjacent-difference vector.

The actual minimum spacing is

    d_actual = min(diff(sort(x)))

and the default binding criterion is

    abs(d_actual-d_min) <= 1e-4.

The complete vectors are additionally printed line-by-line so MATLAB table-width
formatting cannot hide or truncate M=32 constellation entries.

No global-optimum claim is made.
