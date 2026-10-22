function paths = plot_commit_o_m8_snr30_fig7(T,outDir)
%PLOT_COMMIT_N_M16_FIG7  Identical Commit M design for M=8 SNR=30 dB.
% No error bars, scatter points, std annotations, or numeric recomputation.
    names=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    assert(istable(T)&&height(T)==4 && ...
        all(ismember({'Method','ValidatedAMI'},T.Properties.VariableNames)) && ...
        isequal(string(T.Method(:)),names), ...
        'CommitO:PlotSchema','Expected frozen baseline + GPS + SA, in original order.');
    vals=double(T.ValidatedAMI(:));
    assert(all(isfinite(vals))&&all(vals>=0)&&all(vals<=log2(8)+1e-6), ...
        'CommitO:PlotValues','Invalid validated AMI.');
    if ~isfolder(outDir),mkdir(outDir);end
    fig=figure('Color','w','Name','Commit O: Figure 7 M16', ...
        'Units','inches','Position',[1 1 3.55 2.75]);
    ax=axes(fig);hold(ax,'on');box(ax,'on');grid(ax,'on');
    b=bar(ax,1:4,vals,0.62,'FaceColor','flat','EdgeColor',[.3 .3 .3]);
    b.CData=[.78 .78 .78;.61 .61 .61;.43 .43 .43;.22 .22 .22];
    % Intentionally NO errorbar(): statistics of the 12 are tabulated separately.
    ylim(ax,[0 max(vals)*1.15]);
    xlim(ax,[0.5 4.5]);xticks(ax,1:4);
    xticklabels(ax,{'Uniform PAM','Zhou GS','Pattern Search','SA'});
    xtickangle(ax,12);
    xlabel(ax,'Constellation-design method');
    ylabel(ax,'Validated AMI [bits/symbol]');
    title(ax,'M = 8, SNR = 30 dB, \sigma_R^2 = 0.3', ...
        'Interpreter','tex','FontWeight','normal');
    for k=1:4
        text(ax,k,vals(k)+0.018,sprintf('%.4f',vals(k)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',7.2,'Color','k');
    end
    ax.FontName='Times New Roman';ax.FontSize=8;ax.LineWidth=.7;
    ax.GridAlpha=.23;ax.XGrid='off';
    try,ax.Toolbar.Visible='off';catch,end
    drawnow;
    stem='Fig7_CommitO_M8_SNR30_maxAMI_sigmaR2_0p3';
    paths=struct('pdf',fullfile(outDir,[stem '.pdf']), ...
        'png',fullfile(outDir,[stem '.png']), ...
        'fig',fullfile(outDir,[stem '.fig']));
    savefig(fig,paths.fig);
    exportgraphics(fig,paths.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,paths.pdf,'ContentType','vector','BackgroundColor','white');
    catch
        exportgraphics(fig,paths.pdf,'BackgroundColor','white');
    end
    close(fig);
end
