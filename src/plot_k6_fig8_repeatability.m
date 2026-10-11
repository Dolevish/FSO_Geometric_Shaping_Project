function paths = plot_k6_fig8_repeatability(T,statistics,varargin)
%PLOT_K6_FIG8_REPEATABILITY  Twenty-column validated SA restart bar chart.
%
% Each column is ONE independently validated SA trajectory (not the
% winner across multiple restarts). Dashed line: sample mean across all 20.
% Report sample variance (denominator 19), std and observed range.
% Main y-axis begins at zero. A zoomed inset appears only if the observed
% range is too small to be distinguished reliably on the zero-baseline plot.
%
% No SA, fast AMI evaluation or independent validation occurs here.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'CloseFigure',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    if ~istable(T)||height(T)~=20 || ...
            ~all(ismember({'Restart','ValidatedAMI'}, ...
                T.Properties.VariableNames)) || ...
            ~isequal(T.Restart(:),(1:20).')
        error('plot_k6_fig8_repeatability:Table', ...
            'Expected 20 validated AMI rows ordered by restart ID 1:20.');
    end
    y=double(T.ValidatedAMI(:));
    if any(~isfinite(y))||any(y<0)
        error('plot_k6_fig8_repeatability:AMI','Non-finite/negative AMI.');
    end
    req={'MeanAMI','SampleVarianceAMI','SampleStdAMI', ...
        'MinAMI','MaxAMI','RangeAMI'};
    if ~isstruct(statistics)||~all(isfield(statistics,req))
        error('plot_k6_fig8_repeatability:Stats','Missing sample statistics.');
    end
    if abs(mean(y)-statistics.MeanAMI)>1e-10 || ...
            abs(var(y,0)-statistics.SampleVarianceAMI)>1e-10 || ...
            abs(std(y,0)-statistics.SampleStdAMI)>1e-10 || ...
            abs((max(y)-min(y))-statistics.RangeAMI)>1e-10
        error('plot_k6_fig8_repeatability:StatsMismatch', ...
            'The 20 bars and the statistical summary do not agree.');
    end

    if strlength(string(o.OutputDirectory))==0
        error('plot_k6_fig8_repeatability:OutputDirectory', ...
            'Provide an explicit output directory under results/ieee_revision2.');
    end
    outDir=char(string(o.OutputDirectory));
    root=char(java.io.File(fso_result_utils.revision2_results_root()).getCanonicalPath());
    target=char(java.io.File(outDir).getCanonicalPath());
    if ~(strcmp(target,root)||startsWith(target,[root filesep]))
        error('plot_k6_fig8_repeatability:ResultRoot', ...
            'Figure output must remain under results/ieee_revision2.');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('plot_k6_fig8_repeatability:mkdir','%s',msg);end
    end

    fig=figure('Color','w','Name','Commit K6 - Fig. 8 SA repeatability', ...
        'Units','inches','Position',[1 1 5.55 3.2]);
    ax=axes(fig,'Position',[.105 .20 .86 .71]);
    hold(ax,'on');box(ax,'on');grid(ax,'on');

    palette=repmat([.36 .52 .66],20,1);
    [~,iBest]=max(y);
    palette(iBest,:)=[.11 .29 .43];

    b=bar(ax,1:20,y,0.78,'FaceColor','flat', ...
        'EdgeColor',[.16 .28 .40],'LineWidth',.35);
    b.CData=palette;
    hMean=yline(ax,statistics.MeanAMI,'--', ...
        'LineWidth',1.3,'Color',[.73 .28 .14]);
    % Keep the dashed mean clearly visible over the bars.
    try,uistack(hMean,'top');catch,end

    xlim(ax,[.35 20.65]);
    ymax=max(y);
    ylim(ax,[0 max(1.05*ymax,1e-6)]);
    xticks(ax,1:20);
    xticklabels(ax,string(1:20));
    xlabel(ax,'Independent SA restart index');
    ylabel(ax,'Validated AMI [bits/symbol]');
    title(ax,'M = 8, SNR = 20 dB, \sigma_R^2 = 0.1', ...
        'Interpreter','tex','FontWeight','normal','FontSize',10);
    ax.FontName='Times New Roman';
    ax.FontSize=8;
    ax.LineWidth=.8;
    ax.GridColor=[.80 .80 .80];
    ax.GridAlpha=.45;
    ax.XGrid='off';
    try,ax.Toolbar.Visible='off';catch,end

    text(ax,.015,.98,sprintf( ...
        'Mean: %.6f   Std: %.3g   Variance: %.3g   Range: %.3g', ...
        statistics.MeanAMI,statistics.SampleStdAMI, ...
        statistics.SampleVarianceAMI,statistics.RangeAMI), ...
        'Units','normalized','FontName','Times New Roman', ...
        'FontSize',8.0,'VerticalAlignment','top', ...
        'HorizontalAlignment','left', ...
        'BackgroundColor','white','Margin',2,'EdgeColor',[.88 .88 .88]);

    % Optional zoom inset if zero-baseline bars cannot convey the spread.
    if statistics.RangeAMI<0.075 && statistics.RangeAMI>0
        za=axes(fig,'Position',[.58 .48 .35 .25], ...
            'Color','white','Box','on');
        hold(za,'on');grid(za,'on');
        bar(za,1:20,y,0.72,'FaceColor',[.40 .59 .74], ...
            'EdgeColor','none');
        yline(za,statistics.MeanAMI,'--','Color',[.73 .28 .14], ...
            'LineWidth',.8);
        pad=max(.003,statistics.RangeAMI*.18);
        ylim(za,[statistics.MinAMI-pad statistics.MaxAMI+pad]);
        xlim(za,[.35 20.65]);
        xticks(za,[1 5 10 15 20]);
        title(za,'Zoomed variation','FontWeight','normal','FontSize',7.5);
        za.FontName='Times New Roman';
        za.FontSize=6.6;
        za.LineWidth=.7;
        za.GridAlpha=.25;
        za.XGrid='off';
        try,za.Toolbar.Visible='off';catch,end
    end

    drawnow;
    stem='Fig8_SA_repeatability';
    paths=struct();
    paths.pdf=fullfile(outDir,[stem '.pdf']);
    paths.png=fullfile(outDir,[stem '.png']);
    paths.fig=fullfile(outDir,[stem '.fig']);
    savefig(fig,paths.fig);
    exportgraphics(fig,paths.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,paths.pdf,'ContentType','vector', ...
            'BackgroundColor','white');
    catch
        exportgraphics(fig,paths.pdf,'BackgroundColor','white');
    end
    fprintf('K6 Fig.8 saved:\n  %s\n  %s\n  %s\n', ...
        paths.pdf,paths.png,paths.fig);
    if o.CloseFigure,close(fig);end
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
