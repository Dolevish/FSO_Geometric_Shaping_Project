# Commit K — canonical minimum-spacing sensitivity rerun

Commit K resets the minimum-spacing analysis using the final numerical system rather
than the historical experimental scripts under `min_gap_test`.

## Scientific question

How sensitive is the validated optimized AMI to the numerical regularization
constraint d_min when all other channel, evaluator and SA settings are frozen?

## Frozen grid

- M = {16, 32}
- SNR = 20 dB
- sigma_R^2 = 0.2
- d_min = {0, 0.005, 0.01, 0.02, 0.05}
- 2 independent replicates
- 8 restarts per replicate
- 8000 SA iterations per restart for both M values

Totals:

- 10 physical sensitivity points
- 20 replicate optimization cases
- 160 independently validated restart tasks
- 1,280,000 SA iterations
- 1,280,160 nominal fast-objective calls including each restart's initial score

The M=16 budget is intentionally raised from the I.11.1 standard 4000 iterations to
8000 so M=16 and M=32 are compared under the same conservative budget.

## Frozen evaluator and SA system

Commit K inherits the I.11.1 production configuration:

- log-h fast fading evaluator
- candidate-adaptive nonuniform y grid
- YGridScale = 1.1
- log-h span = 8 sigma
- y block size = 512
- case-adaptive I.10.3 dt policy
- independent adaptive-integration AMI validator
- T0 = 0.4
- Tf = 1e-3
- base proposal std = 0.15
- 50 iterations per temperature block
- the same adaptive proposal-step policy
- P_avg = 1
- pinZero = true
- no artificial x_max constraint in the scientific feasible set

At SNR=20 dB and sigma_R^2=0.2, the I.10.3 rule selects dt=0.01 for both M=16
and M=32.

## Paired randomization across d_min

For a fixed (M, replicate), all five d_min values use the same case master seed.
Restart r also uses the same Threefry substream r. The projector changes with d_min,
as it must, but the random construction is paired so a spacing value is not given a
different set of random starts by accident.

## Validation and winner selection

Every restart winner is independently validated. For each
(M, d_min, replicate), the final winner is selected by validated AMI across all eight
restarts, not by fast AMI.

The output reports whether the fast winner and validated winner differ.

## Parallel scheduler

The run uses an 8-process dynamic parfeval/fetchNext queue. One restart is one task.
M=32 work is dispatched first, with smaller d_min values prioritized within M.
Whenever a worker completes a restart, the next pending task is submitted immediately.
The run records queue occupancy and worker utilization.

## Checkpoint/resume

Each restart is atomically checkpointed with the complete Commit-K scientific
signature and task identity. A resumed run reuses only exact matches.

## Results-root policy

Commit K introduces the clean revision-2 root:

    results/ieee_revision2

Historical `results/ieee_revision` artifacts are not moved or overwritten. Commit K
and subsequent clean reruns should use `fso_result_utils.revision2_results_root()`.

Commit-K outputs are stored under:

    results/ieee_revision2/commit_k_min_gap/K_<UTC timestamp>/

## Final outputs

- K_min_gap_results.mat
- K_min_gap_restart_results.csv
- K_min_gap_replicate_summary.csv
- K_min_gap_physical_summary.csv
- K_scheduler_trace.csv
- K_run_manifest.mat
- checkpoints/
- Fig_K_min_gap_sensitivity.pdf
- Fig_K_min_gap_sensitivity.png
- Fig_K_min_gap_sensitivity.fig

The figure contains one set of axes with two curves, M=16 and M=32. Each point is
the mean validated AMI of the two replicate winners; error bars are the sample
standard deviation across replicates.

No global-optimum claim is made.
