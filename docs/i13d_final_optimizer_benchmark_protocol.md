# I.13D Final Optimizer Benchmark — Frozen Protocol

This benchmark is frozen before the final Pattern Search results are observed.

## Scientific question

Are the validated AMI values reported by the frozen I.11.1 Simulated Annealing
production materially dependent on the choice of optimizer?

## Compared optimizers

- SA: frozen I.11.1 production results; they are not rerun or retuned.
- Pattern Search: MATLAB Global Optimization Toolbox, direct linear constraints,
  same fast AMI objective and independent validator.

No global-optimum claim is made for either optimizer.

## Benchmark grid

- M = {4, 8, 16, 32}
- SNR = {20, 30} dB
- sigma_R^2 = {0, 0.3}
- replicates = {1, 2}

Total: 16 physical points, 32 optimization cases.

The grid intentionally contains all constellation sizes, no turbulence vs the strongest
tested turbulence, medium/high SNR, and both independent production replicates.

## Fairness policy

For each case Pattern Search uses:

- the exact same raw start for restart k as frozen SA restart k;
- the same Threefry master seed/substream construction;
- the same average-power, non-negativity, zero-anchor and d_min=0.01 constraints;
- the same fast AMI evaluator;
- the same independent validator;
- the same maximum objective-evaluation budget:
  - M=4,8,16: 4001 evaluations/start;
  - M=32: 8001 evaluations/start.
- no separate pre-score of the Pattern Search initial point is performed;
  its initial evaluation is part of MaxFunctionEvaluations, so no hidden
  extra AMI call is added for logging.

UseCompletePoll is frozen to false so Pattern Search is not intentionally forced to
finish an entire poll after an improvement when the comparison is evaluation-budget
limited. Actual search evaluations are recorded and any overrun is reported.

Because Pattern Search enforces linear constraints only up to numerical solver
tolerance, its returned point may miss an active d_min boundary by around 1e-6.
After termination only, I.13D applies the project's canonical constellation projector
once. The repair must satisfy ||dx||_inf <= 1e-5 or the task fails. This repair does
not feed back into the optimizer. One post-hoc fast AMI audit evaluation and the
independent validator then score the repaired feasible point; the audit evaluation is
reported separately and is not counted as a search evaluation.

## Pre-declared equivalence threshold

A physical point is classified as optimizer-equivalent when

    |mean_rep I_SA - mean_rep I_PS| <= 1e-3 bit/symbol.

Otherwise the sign of I_SA-I_PS determines which optimizer achieved the higher
validated AMI.

## Parallel scheduler

The final run does not use parfor.

Task granularity is one atomic restart:

    (M, SNR, sigma_R^2, replicate, restart).

Eight process workers are maintained by a dynamic parfeval queue. Exactly one
restart future is assigned per worker. On every completion, fetchNext is followed
immediately by dispatch of the next pending restart. Therefore there is no persistent
parfor chunk ownership and no worker can hold a private chain of future restarts.

The queue is longest-first and gives priority to refined-dt turbulent tasks, then
other turbulent tasks, then high M/SNR work.

Pattern Search internal parallelism is disabled; one worker means one restart.

The run records both:

1. queue occupancy — whether all 8 futures remained in flight while pending work
   existed; and
2. measured worker utilization
       sum(task wall seconds)/(workers * scheduler wall seconds).

The unavoidable final tail begins only after fewer than 8 restart tasks remain.

## Checkpointing

Every successful restart is atomically checkpointed on the worker. Resume validates
a complete scientific signature before reusing any checkpoint. The final reduction
is performed only after all 224 restart results are available.

## Saved outputs

- i13d_optimizer_benchmark.mat
- i13d_restart_tasks.csv
- i13d_case_summary.csv
- i13d_physical_summary.csv
- i13d_scheduler_trace.csv
- Fig_I13D_SA_vs_PatternSearch.{pdf,png,fig}
- Fig_I13D_four_method_comparison.{pdf,png,fig} when I.13B Zhou data are available
- Fig_I13D_restart_variability.{pdf,png,fig}

The four-method figure compares Uniform PAM, Zhou-2025 literature GS, Pattern Search
channel-adaptive GS, and frozen-SA channel-adaptive GS on the benchmark grid.
