# Commit O — M8/SNR30 high-SNR optimizer and GPS mesh sensitivity

## Scientific design
- Operating point: **M=8**, **SNR=30 dB**, **sigma_R^2=0.3**, **Pavg=1**, **x1=0**, adjacent gap **0.01**, equiprobable symbols, no CSI.
- Exactly **12 distinct feasible initial constellations**, created via one deterministic Threefry master seed (uniform PAM first; 11 Gaussian-perturbed starts with scale 0.2). **Same paired start per restart index** across all 3 arms; 36 tasks total.
- **SA**: frozen I11.1 settings, 4000 iterations/restart (4001 evaluations including initialization), same T0/Tf/step adaptation.
- **GPS (main)**: complete poll, classical GPSPositiveBasis2N as provided by original MATLAB patternsearch, **MeshTolerance=1e-7**, max 4001 objective evaluations, StepTolerance=1e-8 and FunctionTolerance=1e-8.
- **GPS_CONTROL**: identical to main GPS except **MeshTolerance=1e-6**; same 12 starts, same budget, same validation. This directly isolates mesh tolerance at M8/SNR30.
- Fast objective: original Log-Normal no-CSI Fast AMI with **loghDt=0.01** (case-adaptive policy selects base spacing for M8/SNR30), +/-8 sigma fading span, adaptive y-grid (scale 1.1). Each final constellation independently validated.
- Construct **Uniform PAM** and **Zhou et al. GS** for M8, with Zhou regularization **ExtraBits=2**, and freshly compute both their validated AMIs at SNR=30 dB (no copied SNR20 values).
- Main **Figure 7** shows four grayscale bars Uniform PAM, Zhou GS, tight GPS, SA, in same K5.2/Commit M visual style, **without error bars**. For the two optimizers use **MAX over 12 validated restart AMIs**. The 12-run sample mean, variance and standard deviation are saved separately for all three arms.
- Save the paired GPS-control differences per restart (validated AMI, evaluations, exitflag and final mesh size) and their summary statistics.
- Compare all four bars and optimizer statistics across Commit M (M8/SNR20/mesh1e-6), Commit N (M16/SNR20/mesh1e-6), Commit O (M8/SNR30/mesh1e-7). **The cross-commit variation does not isolate SNR because mesh tolerance also changes; use within-O GPS vs GPS_CONTROL for the mesh effect.**
- **Every worker** forwards its actual **START** and **DONE** through parallel.pool.DataQueue, with **elapsed time in minutes and seconds** such as `04m 32.16s`. Baseline validation prints the same duration format.
- Resume is protected by an exact scientific signature and task-specific checkpoints; never overwrite M/N historical data.

## Run on local MATLAB (branch ieee-revision)
```matlab
addpath(genpath(fullfile(pwd,'src')));
test_commit_o_m8_snr30_preflight();
plan=run_commit_o_m8_snr30('DryRun',true);
O=run_commit_o_m8_snr30();
```
Optional custom historical paths:
```matlab
O=run_commit_o_m8_snr30('M8ResultPath','/exact/path/CommitM_results.mat', ...
    'M16ResultPath','/exact/path/CommitN_M16_results.mat');
```
Resume incomplete run with `'ResumeDirectory','/exact/path/N...'` replaced by the actual **O_...** run directory. Optional `'UseParallel',false` for debugging; default 8 process workers.

## Output after numerical execution
```text
results/ieee_revision2/commit_o_m8_snr30/O_<timestamp>/
  CommitO_all_36_restarts.csv
  CommitO_optimizer_stats_12_restarts.csv
  Fig7_CommitO_four_methods.csv
  Fig7_CommitO_M8_SNR30_maxAMI_sigmaR2_0p3.pdf
  Fig7_CommitO_M8_SNR30_maxAMI_sigmaR2_0p3.png
  Fig7_CommitO_M8_SNR30_maxAMI_sigmaR2_0p3.fig
  CommitO_GPS_mesh_12_pairs.csv
  CommitO_GPS_mesh_summary.csv
  CommitO_M_N_O_four_methods.csv
  CommitO_M_N_O_optimizer_stats.csv
  CommitO_results.mat
  checkpoints/*.mat
  native_gps/*.mat
```
No new MATLAB optimization results are claimed by this code-only commit. All measured figures require running this script locally.