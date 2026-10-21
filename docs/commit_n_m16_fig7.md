# Commit N — M16 paired optimizer comparison

Uses the same scientific configuration as Commit M except M=16: M-PAM positive, SNR 20dB, sigma_R^2=0.3, average optical intensity 1, d_min=0.01, x1=0. 12 unique paired initial constellations (Uniform PAM first and 11 Threefry Gaussian perturbations with standard deviation 0.2); same raw/projected starting point for SA and GPS of the same restart index.

SA: original 4000 iterations per restart, original temperature and step parameters. GPS: classic GPSPositiveBasis2N, UseCompletePoll=true, same I13D mesh/step/function tolerances, 4001 max evaluations per restart. Objective is the same fast log-h AMI (dt=0.01, +/-8 sigma_t fading span, adaptive y-grid scale 1.1); all candidates are independently validated post-optimization.

Each of 24 worker tasks emits START and DONE with elapsed seconds and validated AMI. Messages from multiple workers can interleave. Progress is forwarded to the MATLAB client via parallel.pool.DataQueue.

Both 16-PAM Uniform and Zhou-2025 GS baselines are constructed afresh for M16 and independently evaluated with the same AMI validator. Zhou uses ExtraBits=2 with Regularize=true and is not re-optimized.

Figure 7 for M16 uses the same four grayscale bars and formatting as Commit M (no error bars); optimizer bars show the maximum of 12 independently validated AMIs. Descriptive means, sample variances and sample standard deviations of all 12 are retained in the CSV. The M8 comparison CSV is made from the completed CommitM_results.mat (automatic lookup, or supply M8ResultPath). Both raw AMI and normalized I/log2(M) are reported; cross-M absolute gains do not imply better equal-rate performance.

## Execute from local repository root (branch ieee-revision)

    addpath(genpath(fullfile(pwd,'src')));
    test_commit_n_m16_preflight();
    plan = run_commit_n_m16_fig7('DryRun',true);
    N = run_commit_n_m16_fig7();

If M8 result is stored elsewhere:

    N = run_commit_n_m16_fig7('M8ResultPath','/exact/path/CommitM_results.mat');

Resume an interrupted run:

    N = run_commit_n_m16_fig7('ResumeDirectory','/exact/results/ieee_revision2/commit_n_m16_fig7/N_timestamp');

By default the job uses eight process workers. Serial debug is available with 'UseParallel',false.

## Output, generated after a successful MATLAB run

results/ieee_revision2/commit_n_m16_fig7/N_timestamp/

- CommitN_M16_all_24_restarts.csv
- CommitN_M16_optimizer_stats_12_restarts.csv
- Fig7_CommitN_M16_four_methods.csv
- CommitN_M8_vs_M16_four_methods.csv
- CommitN_M8_vs_M16_optimizer_stats.csv
- Fig7_CommitN_M16_maxAMI_sigmaR2_0p3.pdf/.png/.fig
- CommitN_M16_results.mat
- checkpoints and native_gps task data

Scientific status: source-only commit, not newly measured numerical results. Execution must occur on user's MATLAB installation with Global Optimization and Parallel Computing Toolboxes. No global-optimum claim.