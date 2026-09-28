function bundle = sim_shapingGain_vs_turbulence(M, SNR_dB, varargin)
%SIM_SHAPINGGAIN_VS_TURBULENCE  Validated GS gain versus turbulence.
%
%   bundle = sim_shapingGain_vs_turbulence(M, SNR_dB)
%
% Uses canonical configuration, Commit-C validated multi-start selection,
% deterministic seeds and the same result schema as sim_AMI_vs_SNR.m.
%
% Defaults preserve the original focused experiment: M=16, SNR=20 dB.

    if nargin < 1 || isempty(M), M = 16; end
    if nargin < 2 || isempty(SNR_dB), SNR_dB = 20; end

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'M', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=2 && mod(v,1)==0);
    addRequired(p, 'SNR_dB', @(v) isnumeric(v) && isscalar(v) && isfinite(v));
    addParameter(p, 'P_avg', 1, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p, 'TurbulenceVec', [0 0.1 0.2 0.3], @(v) isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>=0));
    addParameter(p, 'ghN_h', 40, @(v) isnumeric(v) && isscalar(v) && v>=2 && mod(v,1)==0);
    addParameter(p, 'saMaxIter', 50000, @(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p, 'saNStarts', 6, @(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p, 'saUseParallel', true, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'minGap', 0.05, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p, 'pinZero', true, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'baseSeed', 20260928, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p, 'historyEvery', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0));
    addParameter(p, 'ComputeBER', false, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'ResultsRoot', fso_result_utils.default_results_root(), @(v) ischar(v) || isstring(v));
    parse(p, M, SNR_dB, varargin{:});
    o = p.Results;

    M          = double(M);
    SNR_dB     = double(SNR_dB);
    P_avg      = double(o.P_avg);
    sigVec     = double(o.TurbulenceVec(:).');
    nSig       = numel(sigVec);
    xMaxBound  = fso_result_utils.feasible_xmax_bound(M, P_avg, o.minGap, o.pinZero);

    resultsRoot = char(o.ResultsRoot);
    outDir = fullfile(resultsRoot, 'turbulence_sweep', sprintf('M%d_SNR%g', M, SNR_dB));
    caseDir = fullfile(outDir, 'cases');
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end

    runMeta = fso_result_utils.run_metadata(mfilename);
    caseResults = cell(1,nSig);
    caseFiles   = cell(1,nSig);

    fprintf('\n============================================================\n');
    fprintf('IEEE revision turbulence sweep\n');
    fprintf('M=%d | SNR=%.2f dB | Pavg=%.4g | minGap=%.4g\n', M, SNR_dB, P_avg, o.minGap);
    fprintf('SA: %d restarts x %d iterations | restart parallelism=%d\n', ...
        o.saNStarts, o.saMaxIter, o.saUseParallel);
    fprintf('Results: %s\n', outDir);
    fprintf('============================================================\n');

    for idx = 1:nSig
        sigma_X_sq = sigVec(idx);
        seed = fso_result_utils.case_seed(o.baseSeed, M, SNR_dB, sigma_X_sq);

        cfg = build_fso_config(M, P_avg, SNR_dB, sigma_X_sq, ...
            'ghN_h', o.ghN_h, ...
            'xMaxBound', xMaxBound, ...
            'saMaxIter', o.saMaxIter, ...
            'saNStarts', o.saNStarts, ...
            'saUseParallel', o.saUseParallel, ...
            'minGap', o.minGap, ...
            'pinZero', o.pinZero, ...
            'seedInit', seed, ...
            'logEvery', 5000, ...
            'historyEvery', o.historyEvery);

        xPAM = define_constellation(M, P_avg, true, "mean");
        AMI_functions.assert_constellation_feasible(xPAM, cfg, 1e-10);

        amiPAMFast = cfg.AMI_Evaluator(xPAM);
        tBaseVal = tic;
        amiPAMVal = cfg.AMI_Validator(xPAM);
        baseValRuntime = toc(tBaseVal);

        fprintf('\n[%d/%d] sigma_X^2=%.4f | baseline validated AMI=%.6f\n', ...
            idx, nSig, sigma_X_sq, amiPAMVal);

        [saOut, starts] = sa_multistart(cfg, xPAM, ...
            sprintf('[M%d|SNR%.1f|sig%.3f]',M,SNR_dB,sigma_X_sq));

        xGS       = saOut.bestXValidated(:);
        amiGSVal  = saOut.bestMIValidated;
        amiGSFast = saOut.validatedWinnerFastMI;

        berPAM = NaN;
        berGS  = NaN;
        if o.ComputeBER
            berPAM = BER_functions.calculate_BER_noCSI_ML(xPAM(:), cfg);
            berGS  = BER_functions.calculate_BER_noCSI_ML(xGS(:), cfg);
        end

        saOutSaved = saOut;
        if isfield(saOutSaved,'bestRun'), saOutSaved = rmfield(saOutSaved,'bestRun'); end

        c = struct();
        c.meta = runMeta;
        c.case = struct('M',M, 'P_avg',P_avg, 'SNR_dB',SNR_dB, ...
            'sigma_X_sq',sigma_X_sq, 'seed',seed);
        c.config = fso_result_utils.config_snapshot(cfg);
        c.baseline = struct('x',xPAM(:), 'amiFast',amiPAMFast, ...
            'amiValidated',amiPAMVal, 'ber',berPAM);
        c.starts = starts;
        c.optimizer = saOutSaved;
        c.winner = struct('startIndex',saOut.bestStartValidated, 'x',xGS, ...
            'amiFast',amiGSFast, 'amiValidated',amiGSVal, 'ber',berGS);
        c.selection = struct( ...
            'fastWinnerIndex',saOut.bestStartFast, ...
            'validatedWinnerIndex',saOut.bestStartValidated, ...
            'changed',saOut.selectionChanged, ...
            'validatedAdvantageOverFastWinner',saOut.validatedAdvantageOverFastWinner, ...
            'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap, ...
            'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
        c.metrics = struct('gainBits',amiGSVal-amiPAMVal, ...
            'fastGainBits',amiGSFast-amiPAMFast);
        c.runtime = struct('sa',saOut.totalSARuntime, ...
            'validation',saOut.validationTotalRuntime+baseValRuntime, ...
            'total',saOut.totalRuntime+baseValRuntime);

        caseResults{idx} = c;
        caseFiles{idx} = fullfile(caseDir, ...
            fso_result_utils.case_filename(M,SNR_dB,sigma_X_sq));
        caseResult = c; %#ok<NASGU>
        save(caseFiles{idx}, 'caseResult', '-v7.3');

        fprintf('  validated GS AMI=%.6f | gain=%+.6f | fast winner=%d | validated winner=%d\n', ...
            amiGSVal, c.metrics.gainBits, saOut.bestStartFast, saOut.bestStartValidated);
        fprintf('  constellation:'); fprintf(' %.6f', xGS); fprintf('\n');
    end

    summary = struct();
    summary.amiPAMFast      = nan(1,nSig);
    summary.amiPAMValidated = nan(1,nSig);
    summary.amiGSFast       = nan(1,nSig);
    summary.amiGSValidated  = nan(1,nSig);
    summary.gainBits        = nan(1,nSig);
    summary.selectionChanged = false(1,nSig);
    summary.meanAbsFastValidationGap = nan(1,nSig);
    summary.maxAbsFastValidationGap  = nan(1,nSig);
    summary.saRuntime = nan(1,nSig);
    summary.xGS = cell(1,nSig);
    if o.ComputeBER
        summary.berPAM = nan(1,nSig);
        summary.berGS  = nan(1,nSig);
    end

    for idx = 1:nSig
        c = caseResults{idx};
        summary.amiPAMFast(idx)      = c.baseline.amiFast;
        summary.amiPAMValidated(idx) = c.baseline.amiValidated;
        summary.amiGSFast(idx)       = c.winner.amiFast;
        summary.amiGSValidated(idx)  = c.winner.amiValidated;
        summary.gainBits(idx)        = c.metrics.gainBits;
        summary.selectionChanged(idx)= c.selection.changed;
        summary.meanAbsFastValidationGap(idx) = c.selection.meanAbsFastValidationGap;
        summary.maxAbsFastValidationGap(idx)  = c.selection.maxAbsFastValidationGap;
        summary.saRuntime(idx)       = c.runtime.sa;
        summary.xGS{idx}             = c.winner.x;
        if o.ComputeBER
            summary.berPAM(idx) = c.baseline.ber;
            summary.berGS(idx)  = c.winner.ber;
        end
    end

    bundle = struct();
    bundle.meta       = runMeta;
    bundle.experiment = 'shaping_gain_vs_turbulence';
    bundle.axes       = struct('M',M, 'P_avg',P_avg, 'SNR_dB',SNR_dB, 'sigma_X_sq',sigVec);
    bundle.settings   = struct( ...
        'ghN_h',o.ghN_h, 'saMaxIter',o.saMaxIter, 'saNStarts',o.saNStarts, ...
        'minGap',o.minGap, 'pinZero',o.pinZero, 'baseSeed',o.baseSeed, ...
        'historyEvery',o.historyEvery, 'computeBER',o.ComputeBER, ...
        'xMaxBound',xMaxBound);
    bundle.caseFiles = caseFiles;
    bundle.summary   = summary;

    bundlePath = fullfile(outDir, sprintf('shapingGain_M%d_SNR%g.mat',M,SNR_dB));
    save(bundlePath, 'bundle', '-v7.3');

    fprintf('\nValidated summary\n');
    fprintf(' sigma_X^2 | PAM AMI | GS AMI | gain | selection changed\n');
    for idx = 1:nSig
        fprintf('   %7.4f | %.6f | %.6f | %+.6f | %d\n', ...
            sigVec(idx), summary.amiPAMValidated(idx), ...
            summary.amiGSValidated(idx), summary.gainBits(idx), ...
            summary.selectionChanged(idx));
    end
    fprintf('Saved summary bundle: %s\n', bundlePath);

    plot_bundle(bundle);
end


function plot_bundle(bundle)
    sig = bundle.axes.sigma_X_sq;
    S = bundle.summary;

    figure('Name','Validated AMI vs turbulence','Color','w');
    plot(sig, S.amiPAMValidated, '--o', 'LineWidth',1.5, 'DisplayName','Uniform PAM');
    hold on; grid on; box on;
    plot(sig, S.amiGSValidated, '-s', 'LineWidth',1.8, 'DisplayName','GS (validated winner)');
    xlabel('Normalized intensity variance \sigma_X^2');
    ylabel('AMI [bits/symbol]');
    title(sprintf('M=%d, SNR=%.1f dB',bundle.axes.M,bundle.axes.SNR_dB));
    legend('Location','best');

    figure('Name','Validated shaping gain vs turbulence','Color','w');
    plot(sig, S.gainBits, '-o', 'LineWidth',1.8);
    grid on; box on;
    xlabel('Normalized intensity variance \sigma_X^2');
    ylabel('AMI gain [bits/symbol]');
    title(sprintf('GS gain: M=%d, SNR=%.1f dB',bundle.axes.M,bundle.axes.SNR_dB));
end
