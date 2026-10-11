function F = plot_k7_fig8_repeatability_dual(T,stats,outputDir,varargin)
%PLOT_K7_FIG8_REPEATABILITY_DUAL Two separate IEEE-style Fig.8 bar charts.
%
% Figure A: 20 standalone bars with zero-origin y axis, NO zoom inset.
%           Summary statistics live ABOVE the plotting axes, never over bars.
% Figure B: standalone zoom on the same 20 bars, with explicitly truncated
%           y axis and no inset. Both show the validated-AMI sample mean.
%
% This function is strictly post-processing. No SA or AMI is evaluated.
% Statistical variance/std MUST use 19-denominator sample estimates
% from frozen K6, never population variance.
%
% F = plot_k7_fig8_repeatability_dual(R6.restartTable,R6.statistics,dir);
% F = plot_k7_fig8_repeatability_dual(T,stats,dir,'CloseFigures',true);

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'CloseFigures',false,@(x)islogical(x)&&isscalar(x));
    parse(p,varargin{:}); opt=p.Results;
    validate_inputs(T,stats,outputDir);

    y=double(T.ValidatedAMI(:));
    idx=1:20;
    yMean=stats.MeanAMI;
    [~,bestIdx]=max(y);

    % Identical bar encoding in both panels, highlight best restart only.
    colors=repmat([0.40 0.58 0.72],20,1);
    colors(bestIdx,:)=[0.12 0.34 0.50];
    meanColor=[0.74 0.28 0.15];

    % ------------------------------------------------------------------
    % Figure A: zero-origin bar chart; stats in independent HEADER band.
    % ------------------------------------------------------------------
    mainFig=figure('Color','w', ...
        'Name','K7 Fig. 8 main - zero-origin SA repeatability', ...
        'Units','inches','Position',[1 1 6.4 3.75]);
    mainAx=axes('Parent',mainFig,'Units','normalized', ...
        'Position',[0.105 0.18 0.865 0.59]);
    style_axis(mainAx);
    draw_bars(mainAx,idx,y,colors,meanColor,yMean);
    xlim(mainAx,[0.4 20.6]);
    ylim(mainAx,[0 max(y)*1.047]);
    xticks(mainAx,1:20);
    xlabel(mainAx,'Independent SA restart index');
    ylabel(mainAx,'Validated AMI [bits/symbol]');

    % Two dedicated annotation rows, physically outside the axes.
    draw_header(mainFig,'M = 8, SNR = 20 dB, \sigma_R^2 = 0.1');
    statsText=sprintf(['Mean %.6f    Sample std %.6f    ' ...
        'Sample variance %.3e    Range %.6f'], ...
        stats.MeanAMI,stats.SampleStdAMI, ...
        stats.SampleVarianceAMI,stats.RangeAMI);
    annotation(mainFig,'textbox',[0.105 0.80 0.865 0.078], ...
        'String',statsText,'Interpreter','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FitBoxToText','off', ...
        'FontName','Times New Roman','FontSize',9, ...
        'EdgeColor',[0.83 0.83 0.83], ...
        'BackgroundColor','white');
    mainStem=fullfile(outputDir,'Fig8_SA_repeatability_main');
    mainPaths=export_triplet(mainFig,mainStem);

    % ------------------------------------------------------------------
    % Figure B: independent zoomed bar chart; its own full canvas.
    % ------------------------------------------------------------------
    zoomFig=figure('Color','w', ...
        'Name','K7 Fig. 8 zoom - truncated SA repeatability y axis', ...
        'Units','inches','Position',[1 1 6.4 3.75]);
    zoomAx=axes('Parent',zoomFig,'Units','normalized', ...
        'Position',[0.105 0.18 0.865 0.60]);
    style_axis(zoomAx);
    draw_bars(zoomAx,idx,y,colors,meanColor,yMean);
    xlim(zoomAx,[0.4 20.6]);
    pad=max(0.0008,0.15*stats.RangeAMI);
    low=stats.MinAMI-pad;
    high=stats.MaxAMI+pad;
    if low<=0
        low=max(eps,min(y)*.75);
    end
    if high<=low
        high=low+max(0.001,abs(low)*0.002);
    end
    ylim(zoomAx,[low high]);
    xticks(zoomAx,1:20);
    xlabel(zoomAx,'Independent SA restart index');
    ylabel(zoomAx,'Validated AMI [bits/symbol]');
    draw_header(zoomFig, ...
        'M = 8, SNR = 20 dB, \sigma_R^2 = 0.1  (zoom)');
    annotation(zoomFig,'textbox',[0.105 0.805 0.865 0.075], ...
        'String',sprintf(['Expanded AMI axis (not zero-based).  ' ...
        'Dashed line: mean = %.6f bits/symbol'],yMean), ...
        'Interpreter','none','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FitBoxToText','off', ...
        'FontName','Times New Roman','FontSize',9, ...
        'EdgeColor','none','BackgroundColor','white');
    zoomStem=fullfile(outputDir,'Fig8_SA_repeatability_zoom');
    zoomPaths=export_triplet(zoomFig,zoomStem);

    F=struct();
    F.version='K7-figure-only-v1';
    F.outputDirectory=char(outputDir);
    F.main=mainPaths;
    F.zoom=zoomPaths;
    F.main.figureHandle=mainFig;
    F.main.axesHandle=mainAx;
    F.main.yLimits=ylim(mainAx);
    F.zoom.figureHandle=zoomFig;
    F.zoom.axesHandle=zoomAx;
    F.zoom.yLimits=ylim(zoomAx);
    F.samples=20;
    F.mean=yMean;
    F.sampleStd=stats.SampleStdAMI;
    F.sampleVariance=stats.SampleVarianceAMI;
    F.range=stats.RangeAMI;

    fprintf('\nK7 Fig. 8 MAIN (zero origin, no inset):\n');
    fprintf('  %s\n  %s\n  %s\n',F.main.pdf,F.main.png,F.main.fig);
    fprintf('K7 Fig. 8 ZOOM (separate figure):\n');
    fprintf('  %s\n  %s\n  %s\n',F.zoom.pdf,F.zoom.png,F.zoom.fig);

    if opt.CloseFigures
        close(mainFig);
        close(zoomFig);
        F.main=rmfield(F.main,{'figureHandle','axesHandle'});
        F.zoom=rmfield(F.zoom,{'figureHandle','axesHandle'});
    end
end


function validate_inputs(T,stats,outDir)
    if ~istable(T)||height(T)~=20 || ...
            ~all(ismember({'Restart','ValidatedAMI'}, ...
            T.Properties.VariableNames)) || ...
            ~isequal(double(T.Restart(:)),(1:20).')
        error('K7:RestartTable', ...
            'Require all 20 K6 trajectories in restart-index order 1:20.');
    end
    y=double(T.ValidatedAMI(:));
    if ~isvector(y)||numel(y)~=20 || ...
            any(~isfinite(y))||any(y<0)
        error('K7:ValidatedAMI','Invalid K6 validated AMI values.');
    end
    req={'N','MeanAMI','SampleVarianceAMI','SampleStdAMI', ...
        'MinAMI','MaxAMI','RangeAMI'};
    if ~isstruct(stats)||~all(isfield(stats,req)) || stats.N~=20
        error('K7:Summary','Missing or invalid K6 sample statistics.');
    end
    vals=[stats.MeanAMI,stats.SampleVarianceAMI,stats.SampleStdAMI, ...
        stats.MinAMI,stats.MaxAMI,stats.RangeAMI];
    if any(~isfinite(vals))||any(vals<0) || ...
            abs(mean(y)-stats.MeanAMI)>1e-10 || ...
            abs(var(y,0)-stats.SampleVarianceAMI)>1e-10 || ...
            abs(std(y,0)-stats.SampleStdAMI)>1e-10 || ...
            abs(min(y)-stats.MinAMI)>1e-10 || ...
            abs(max(y)-stats.MaxAMI)>1e-10 || ...
            abs((max(y)-min(y))-stats.RangeAMI)>1e-10
        error('K7:StatisticsMismatch', ...
            'Measured K6 sample data disagree with stored N-1 statistics.');
    end
    if ~(ischar(outDir)||(isstring(outDir)&&isscalar(outDir))) || ...
            strlength(string(outDir))==0
        error('K7:OutputDir','Output directory is required.');
    end
    root=char(java.io.File( ...
        fso_result_utils.revision2_results_root()).getCanonicalPath());
    target=char(java.io.File(char(outDir)).getCanonicalPath());
    if ~(strcmp(root,target)||startsWith(target,[root filesep]))
        error('K7:OutputRoot', ...
            'Figures must be under results/ieee_revision2.');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('K7:mkdir','%s',msg);end
    end
end


function style_axis(ax)
    hold(ax,'on');box(ax,'on');grid(ax,'on');
    ax.FontName='Times New Roman';
    ax.FontSize=9;
    ax.LineWidth=0.75;
    ax.GridColor=[0.81 0.81 0.81];
    ax.GridAlpha=0.35;
    ax.XGrid='off';
    ax.Color='white';
    try,ax.Toolbar.Visible='off';catch,end
end


function draw_bars(ax,x,y,colors,lineColor,yMean)
    b=bar(ax,x,y,0.78,'FaceColor','flat', ...
        'EdgeColor',[0.19 0.33 0.43],'LineWidth',0.3);
    b.CData=colors;
    m=yline(ax,yMean,'--','Color',lineColor,'LineWidth',1.5);
    try,uistack(m,'top');catch,end
end


function draw_header(fig,heading)
    annotation(fig,'textbox',[0.105 0.925 0.865 0.052], ...
        'String',heading,'Interpreter','tex', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'EdgeColor','none','BackgroundColor','none', ...
        'FontName','Times New Roman','FontSize',12, ...
        'FontWeight','normal','FitBoxToText','off');
end


function paths=export_triplet(fig,stem)
    paths=struct();
    paths.fig=[stem '.fig'];
    paths.png=[stem '.png'];
    paths.pdf=[stem '.pdf'];
    drawnow;
    savefig(fig,paths.fig);
    exportgraphics(fig,paths.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,paths.pdf,'ContentType','vector', ...
            'BackgroundColor','white');
    catch
        exportgraphics(fig,paths.pdf,'BackgroundColor','white');
    end
end
