function fig = plot_AMI_vs_SNR_allM(resultsRoot)
%PLOT_AMI_VS_SNR_ALLM  Plot validated AMI directly from saved result bundles.
%
%   fig = plot_AMI_vs_SNR_allM()
%   fig = plot_AMI_vs_SNR_allM(resultsRoot)
%
% No constellation values are hard-coded and AMI is not recomputed here.
% The script discovers summary bundles produced by sim_AMI_vs_SNR.m under:
%   <resultsRoot>/ami_vs_snr/M*/AMI_vs_SNR_M*.mat
%
% Each subplot corresponds to one M and contains uniform-PAM (dashed) and
% validated-GS (solid) curves for every stored turbulence level.

    if nargin < 1 || isempty(resultsRoot)
        resultsRoot = fso_result_utils.default_results_root();
    end
    resultsRoot = char(resultsRoot);
    amiRoot = fullfile(resultsRoot, 'ami_vs_snr');

    files = dir(fullfile(amiRoot, 'M*', 'AMI_vs_SNR_M*.mat'));
    if isempty(files)
        error('plot_AMI_vs_SNR_allM:NoResults', ...
            ['No AMI-vs-SNR bundles found under %s. Run sim_AMI_vs_SNR(M) ' ...
             'for the desired M values first.'], amiRoot);
    end

    bundles = cell(numel(files),1);
    Mvals = nan(numel(files),1);

    for k = 1:numel(files)
        path = fullfile(files(k).folder, files(k).name);
        D = load(path, 'bundle');
        if ~isfield(D,'bundle') || ~isfield(D.bundle,'experiment') || ...
                ~strcmp(D.bundle.experiment,'AMI_vs_SNR')
            error('plot_AMI_vs_SNR_allM:BadBundle', ...
                'File is not an AMI_vs_SNR bundle: %s', path);
        end
        bundles{k} = D.bundle;
        Mvals(k) = D.bundle.axes.M;
    end

    [Mvals, order] = sort(Mvals);
    bundles = bundles(order);

    % If duplicate bundles exist for one M, retain the first discovered file
    % only after sorting. The normal pipeline writes a single fixed filename.
    [Mvals, uniqueIdx] = unique(Mvals, 'stable');
    bundles = bundles(uniqueIdx);

    nM = numel(Mvals);
    nCols = min(2,nM);
    nRows = ceil(nM/nCols);

    fig = figure('Name','Validated AMI vs SNR - all M','Color','w');
    tl = tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');

    for iM = 1:nM
        b = bundles{iM};
        S = b.summary;
        snr = b.axes.SNR_dB;
        sig = b.axes.sigma_X_sq;

        if ~isequal(size(S.amiPAMValidated), [numel(sig), numel(snr)]) || ...
           ~isequal(size(S.amiGSValidated),  [numel(sig), numel(snr)])
            error('plot_AMI_vs_SNR_allM:DimensionMismatch', ...
                'Summary dimensions do not match axes for M=%d.', b.axes.M);
        end

        ax = nexttile(tl);
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        C = lines(numel(sig));

        for iSig = 1:numel(sig)
            plot(ax, snr, S.amiGSValidated(iSig,:), '-o', ...
                'Color',C(iSig,:), 'LineWidth',1.8, ...
                'DisplayName',sprintf('GS, \\sigma_X^2=%.2f',sig(iSig)));
            plot(ax, snr, S.amiPAMValidated(iSig,:), '--', ...
                'Color',C(iSig,:), 'LineWidth',1.3, ...
                'DisplayName',sprintf('PAM, \\sigma_X^2=%.2f',sig(iSig)));
        end

        yline(ax,log2(b.axes.M),':','HandleVisibility','off');
        xlabel(ax,'SNR [dB]');
        ylabel(ax,'AMI [bits/symbol]');
        title(ax,sprintf('M = %d',b.axes.M));
        xlim(ax,[min(snr) max(snr)]);
        ylim(ax,[0 log2(b.axes.M)*1.05]);
        legend(ax,'Location','best','FontSize',8);
    end

    title(tl,'Validated geometric shaping for No-CSI log-normal FSO');

    outPng = fullfile(amiRoot,'AMI_vs_SNR_allM.png');
    outFig = fullfile(amiRoot,'AMI_vs_SNR_allM.fig');
    try
        exportgraphics(fig,outPng,'Resolution',200);
        savefig(fig,outFig);
        fprintf('Saved figure: %s\n',outPng);
    catch ME
        warning('plot_AMI_vs_SNR_allM:SaveFailed', ...
            'Figure was created but could not be saved: %s',ME.message);
    end
end
