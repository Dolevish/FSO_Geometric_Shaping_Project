# Commit M — final SA/GPS comparison for manuscript Figure 7

## Scientific changes from K5.2
- **Operating point remains frozen:** M=8, SNR=20 dB, sigma_R^2=0.3, Pavg=1, first level pinned to zero, minimum gap 0.01.
- **GPS:** MATLAB classic \`GPSPositiveBasis2N\` (confirmed separately on the project MATLAB R2025b installation), \`UseCompletePoll=true\`; all other I13D GPS options unchanged (mesh/step/function tolerances, no nested parallelism, 4001 maximum objective evaluations).
- **SA:** exact frozen I.11.1 algorithm/temperature/proposal/step-adaptation settings, 4000 iterations per restart, with the initial evaluation included for an intended 4001 objective evaluations.
- **Starts:** exactly 12 distinct *projected* starts from one master seed. Start 1 is Uniform PAM; starts 2–12 are baseline + independent Threefry-substream Gaussian perturbations of scale 0.2. **Each corresponding SA/GPS pair has exactly the same raw and feasible initial constellation**. No duplicate baseline from two old replicates.
- **Scoring:** independently validate the best candidate of EVERY restart. Primary reported GPS/SA AMI is **max across the 12 validated restart results**; save their mean, sample variance, sample standard deviation, minimum, winner restart, and actual evaluation counts.
- **Figure 7:** four grayscale columns with *exactly the same K5.2 layout and aesthetics*. Bars 1–2 are frozen *Uniform PAM* and *Zhou GS* results from the previous validated K5.2 four-method CSV; bars 3–4 show the newly measured GPS/SA maxima. **No error bars** on the figure.
- **Traceability:** frozen reference values are loaded from an existing CSV and audited; no Zhou or uniform optimization is repeated. Outputs go to \`results/ieee_revision2/commit_m_final_fig7/M_<timestamp>/\`, with an optional checkpoint resume directory. Original K5.2 results and Figure 7 are never overwritten.

## MATLAB execution (repo branch \`ieee-revision\`)
\`\`\`matlab
addpath(genpath(fullfile(pwd,'src')));
test_commit_m_preflight();
plan = run_commit_m_final_fig7('DryRun',true);
M = run_commit_m_final_fig7();
\`\`\`
If the frozen K5.2 four-method CSV is not found automatically, specify its **exact** path:
\`\`\`matlab
M = run_commit_m_final_fig7('ReferenceCSVPath', ...
    '/exact/path/Fig7_K5p2_four_methods.csv');
\`\`\`
If the MATLAB task is interrupted:
\`\`\`matlab
M = run_commit_m_final_fig7('ResumeDirectory', ...
    '/exact/path/results/ieee_revision2/commit_m_final_fig7/M_<timestamp>');
\`\`\`
A serial diagnostic is available with \`'UseParallel',false\`. Default runs a dynamic queue of 24 tasks over eight process workers, one full independent optimizer restart per task.

## Generated artifacts, *after* MATLAB finishes
\`\`\`text
results/ieee_revision2/commit_m_final_fig7/M_<timestamp>/
  CommitM_all_24_restarts.csv
  CommitM_optimizer_stats_12_restarts.csv
  Fig7_CommitM_four_methods.csv
  Fig7_CommitM_maxAMI_sigmaR2_0p3.pdf
  Fig7_CommitM_maxAMI_sigmaR2_0p3.png
  Fig7_CommitM_maxAMI_sigmaR2_0p3.fig
  CommitM_results.mat
  checkpoints/SA_restart_01.mat ... GPS_restart_12.mat
  native_gps/CommitM_GPS_native_r01.mat ... r12.mat
\`\`\`
\`CommitM_results.mat\` saves the full restart structures and inputs; the CSVs store numeric values and constellation vectors. GPS \`Exitflag\` and \`FunctionEvaluations\` reveal whether it stopped on mesh tolerance or exhausted the nominal cap.

## Interpretation
The means/stds describe the **12 restart outcomes**, not uncertainty of the maximum or a confidence interval. The graph plots the maximum, so it has no symmetric error bars. No global-optimum or statistically significant SA superiority claim follows from this one operating point.

## Reproducibility boundary
This commit contains source, tests and instructions. It does **not** claim new MATLAB outcomes, because the expensive 24-restart optimization must actually be executed in MATLAB R2025b with Global Optimization and Parallel Computing Toolboxes on the user's machine.
