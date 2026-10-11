function paths = plot_fig7_four_methods(T,sigmaR2,outDir,varargin)
%PLOT_FIG7_FOUR_METHODS  Validated four-method AMI comparison for manuscript.
%
% Exactly four methods: Uniform PAM, Zhou-2025 GS, Pattern Search, SA.
% Input is measured/validated data, never optimized or recomputed here.
% Error bars for PS and SA are sample std across the 2 replicate winners.
% Uniform PAM and Zhou are deterministic and have no replicate error bars.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'CloseFigure',false,@(x)islogical(x)&&isscalar(x));
    parse(p,varargin{:}); o=p.Results;

    required={'Method','ValidatedAMI','StdAcrossReplicates'};
    expected=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    if ~istable(T)||height(T)~=4 || ...
            ~isempty(setdiff(required,T.Properties.VariableNames)) || ...
            ~isequal(string(T.Method(:)),expected)
        error('plot_fig7_four_methods:Data', ...
            'Require four ordered validated methods with AMI and replicate std.');
    end
    vals=double(T.ValidatedAMI(:));
    s=double(T.StdAcrossReplicates(:));
    if any(~isfinite(vals))||any(vals<0)||any(~isfinite(s))|| ...
            any(s<0)||any(s(1:2)~=0)|| ...
            ~isfinite(sigmaR2)||abs(sigmaR2-0.1)>1e-12
        error('plot_fig7_four_methods:ScientificInputs', ...
            'Invalid scores, deterministic errors, or sigma_R^2 != 0.1.');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('plot_fig7_four_methods:Output','%s',msg);end
    end

    fig=figure('Color','w','Name','Commit K5: Fig. 7', ...
        'Units','inches','Position',[1 1 3.55 2.75]);
    ax=axes(fig);hold(ax,'on');box(ax,'on');grid(ax,'on');
    b=bar(ax,1:4,vals,0.62,'FaceColor','flat','EdgeColor',[.3 .3 .3]);
    b.CData=[.78 .78 .78;.61 .61 .61;.43 .43 .43;.22 .22 .22];
    errorbar(ax,[3 4],vals(3:4),s(3:4), ...
        'LineStyle','none','Color','k','LineWidth',0.8,'CapSize',4);
    ylim(ax,[0 max(vals+s)*1.15]);
    xlim(ax,[0.5 4.5]);xticks(ax,1:4);
    xticklabels(ax,{'Uniform PAM','Zhou GS','Pattern Search','SA'});
    xtickangle(ax,12);
    xlabel(ax,'Constellation-design method');
    ylabel(ax,'Validated AMI [bits/symbol]');
    title(ax,'M = 8, SNR = 20 dB, \sigma_R^2 = 0.1', ...
        'Interpreter','tex','FontWeight','normal');
    for k=1:4
        text(ax,k,vals(k)+s(k)+0.018,sprintf('%.4f',vals(k)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',7.2,'Color','k');
    end
    ax.FontName='Times New Roman';ax.FontSize=8;ax.LineWidth=.7;
    ax.GridAlpha=.23;ax.XGrid='off';
    try,ax.Toolbar.Visible='off';catch,end
    drawnow;

    stem='Fig7_four_methods_M8_20dB_sigmaR2_0p1';
    paths=struct();
    paths.pdf=fullfile(outDir,[stem '.pdf']);
    paths.png=fullfile(outDir,[stem '.png']);
    paths.fig=fullfile(outDir,[stem '.fig']);
    savefig(fig,paths.fig);
    exportgraphics(fig,paths.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,paths.pdf,'ContentType','vector','BackgroundColor','white');
    catch
        exportgraphics(fig,paths.pdf,'BackgroundColor','white');
    end
    if o.CloseFigure,close(fig);end
end
