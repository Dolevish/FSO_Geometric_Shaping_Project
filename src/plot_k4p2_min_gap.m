function paths = plot_k4p2_min_gap(result,varargin)
%PLOT_K4P2_MIN_GAP  Same K4 combined one-axis style, NO ERROR BARS.
%
% Exactly 16 real validated points (9 at M16, 7 at M32), beginning
% d_min=0.001. Colors, symbols, axes dimensions, grid and fonts match
% plot_k4_combined_min_gap. No interpolation/extrapolation into the
% infeasible M32 region.
%
% Optional exact user-supplied original FIG graphical template:
%   plot_k4p2_min_gap(C,'TemplateFigPath','/path/to/old_no_errorbars.fig')
% This updates only the two corresponding line XData/YData; old source
% FIG is never changed. Default uses the K4 style programmatically.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@is_text);
    addParameter(p,'TemplateFigPath','',@is_text);
    addParameter(p,'CloseFigure',false,@is_logical);
    parse(p,varargin{:});o=p.Results;

    if ischar(result)||isstring(result)
        A=load(char(result),'result');
        if ~isfield(A,'result')
            error('K4p2:PlotMAT','MAT must contain result struct.');
        end
        result=A.result;
    end
    if ~isstruct(result)||~isfield(result,'physicalSummary') || ...
            ~isfield(result,'config')||~istable(result.physicalSummary)
        error('K4p2:PlotInput','Expected validated merged result struct.');
    end
    T=sortrows(result.physicalSummary,{'M','MinGap'});
    for m=[16 32]
        x=T.MinGap(T.M==m);
        expected=[.001 .005 .01 .02 .03 .04 .05];
        if m==16,expected=[expected .07 .1];end
        if numel(x)~=numel(expected)||any(abs(x(:)-expected(:))>1e-12)
            error('K4p2:PlotGrid','Unexpected physical spacing grid for M=%d.',m);
        end
    end
    if any(~isfinite(T.MeanValidatedAMI)) || ...
            any(~isfinite(T.StdAcrossReplicates)) || ...
            any(abs(T.MinGap)<1e-12) || ...
            height(T)~=16
        error('K4p2:PlotValues','Invalid or zero-gap values.');
    end
    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=char(string(result.outputDirectory));
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    if strlength(string(o.TemplateFigPath))>0
        if exist(char(string(o.TemplateFigPath)),'file')~=2
            error('K4p2:TemplateMissing','FIG template not found.');
        end
        fig=openfig(char(string(o.TemplateFigPath)),'invisible');
        set(fig,'Color','w','Visible','on');
        ax=findobj(fig,'Type','axes');
        ax=ax(~arrayfun(@(a)strcmpi(char(a.Tag),'legend'),ax));
        if numel(ax)~=1
            close(fig);
            error('K4p2:TemplateAxes','Expected one plot axis in the supplied FIG.');
        end
        ax=ax(1);
        h=findobj(ax,'Type','line');
        if numel(h)<2
            close(fig);error('K4p2:TemplateLines','Template needs 2 line objects.');
        end
        chosen=gobjects(1,2);
        for mi=1:2
            m=[16 32];target=m(mi);
            for j=1:numel(h)
                name=string(h(j).DisplayName);
                if contains(name,sprintf('%d',target)) && ...
                        ~contains(name,'boundary')
                    chosen(mi)=h(j);break;
                end
            end
        end
        if ~all(isgraphics(chosen))
            close(fig);
            error('K4p2:TemplateIdentification', ...
                'Could not identify M=16 and M=32 line handles by DisplayName.');
        end
        for mi=1:2
            m=[16 32];rr=T(T.M==m(mi),:);
            set(chosen(mi),'XData',rr.MinGap(:).','YData',rr.MeanValidatedAMI(:).');
        end
        % Preserve original axes/legend styling and update only ticks,
        % x limits, data-driven y limits and data of the two series.
        apply_axis_limits(ax,T,result.config.P_avg);
    else
        fig=figure('Color','w', ...
            'Name','Commit K4.2 combined minimum-spacing sensitivity', ...
            'Units','inches','Position',[1 1 7.8 4.75]);
        ax=axes(fig);hold(ax,'on');box(ax,'on');grid(ax,'on');
        colors=[0 0.4470 0.7410;0.8500 0.3250 0.0980];
        markers={'o','s'};

        % Same physical-feasibility shading as the original combined K4.
        yv=[T.MeanValidatedAMI-T.StdAcrossReplicates; ...
            T.MeanValidatedAMI+T.StdAcrossReplicates];
        span=max(max(yv)-min(yv),0.1);
        ylimNew=[min(yv)-0.07*span max(yv)+0.10*span];
        d32max=2*result.config.P_avg/31;
        patch(ax,[d32max 0.104 0.104 d32max], ...
            [ylimNew(1) ylimNew(1) ylimNew(2) ylimNew(2)], ...
            [0.78 0.78 0.78],'FaceAlpha',0.15,'EdgeColor','none', ...
            'HandleVisibility','off');
        xline(ax,d32max,':','Color',[0.45 0.45 0.45], ...
            'LineWidth',1,'HandleVisibility','off');

        for mi=1:2
            m=[16 32];
            rr=T(T.M==m(mi),:);
            plot(ax,rr.MinGap,rr.MeanValidatedAMI, ...
                'LineStyle','-','Color',colors(mi,:), ...
                'LineWidth',1.5,'Marker',markers{mi}, ...
                'MarkerSize',6,'MarkerFaceColor','w', ...
                'DisplayName',sprintf('M = %d',m(mi)));
        end
        xlabel(ax,'Minimum spacing d_{min}','Interpreter','tex');
        ylabel(ax,'Validated AMI [bits/symbol]');
        title(ax,'SNR = 20 dB, \sigma_R^2 = 0.2', ...
            'Interpreter','tex','FontWeight','normal');
        legend(ax,'Location','northeast');
        ax.Color='w';ax.XColor='k';ax.YColor='k';
        ax.FontName='Times New Roman';ax.FontSize=10;ax.LineWidth=.8;
        ax.GridColor=[.82 .82 .82];ax.GridAlpha=.5;
        apply_axis_limits(ax,T,result.config.P_avg);
    end
    try,ax.Toolbar.Visible='off';catch,end
    drawnow;

    stem='Fig_K4p2_combined_min_gap_no_errorbars';
    paths=struct('fig',fullfile(outDir,[stem '.fig']), ...
        'pdf',fullfile(outDir,[stem '.pdf']), ...
        'png',fullfile(outDir,[stem '.png']));
    savefig(fig,paths.fig);
    exportgraphics(fig,paths.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,paths.pdf,'ContentType','vector','BackgroundColor','white');
    catch
        exportgraphics(fig,paths.pdf,'BackgroundColor','white');
    end
    fprintf('K4.2 combined no-errorbars plot:\n  %s\n  %s\n  %s\n', ...
        paths.pdf,paths.png,paths.fig);
    if o.CloseFigure,close(fig);end
end

function apply_axis_limits(ax,T,Pavg)
    yy=[T.MeanValidatedAMI-T.StdAcrossReplicates; ...
        T.MeanValidatedAMI+T.StdAcrossReplicates];
    span=max(max(yy)-min(yy),.1);
    ylim(ax,[min(yy)-.07*span max(yy)+.1*span]);
    xlim(ax,[-.003 .104]);
    xticks(ax,[.001 .005 .01 .02 .03 .04 .05 .07 .10]);
    xticklabels(ax,{'0.001','0.005','0.01','0.02','0.03', ...
        '0.04','0.05','0.07','0.10'});
    % The shaded region is only a visual indicator; it does not add
    % fictitious post-feasibility measurements.
    assert(Pavg==1,'K4p2:PlotPower','Unexpected mean optical intensity.');
end
function tf=is_text(x),tf=ischar(x)||(isstring(x)&&isscalar(x));end
function tf=is_logical(x),tf=islogical(x)&&isscalar(x);end
