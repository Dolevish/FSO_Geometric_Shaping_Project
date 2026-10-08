function paths = plot_k4_combined_min_gap(result,varargin)
%PLOT_K4_COMBINED_MIN_GAP  Publication-ready merged spacing-sensitivity plot.
%
% Plots mean validated AMI +/- sample std of the 2 independent replicate
% winners at all 16 feasible physical points (M16:9; M32:7).
% Does not extrapolate to infeasible M32 gaps 0.07 and 0.10.
%
%   paths=plot_k4_combined_min_gap(C4);
%   paths=plot_k4_combined_min_gap('/path/to/K4_combined_results.mat');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'CloseFigure',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    if ischar(result)||isstring(result)
        S=load(char(result),'result');
        if ~isfield(S,'result')
            error('plot_k4_combined_min_gap:MAT','MAT file must contain result.');
        end
        result=S.result;
    end
    if ~isstruct(result)||~isfield(result,'physicalSummary') || ...
            ~istable(result.physicalSummary)
        error('plot_k4_combined_min_gap:Input', ...
            'Expected merged K4 result struct with physicalSummary.');
    end

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=char(string(result.outputDirectory));
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    P=sortrows(result.physicalSummary,{'M','MinGap'});
    if height(P)~=16||nnz(P.M==16)~=9||nnz(P.M==32)~=7
        error('plot_k4_combined_min_gap:Grid','Expected 9 M16 + 7 M32 points.');
    end

    fig=figure('Color','w','Name','K4 combined minimum-spacing sensitivity', ...
        'Units','inches','Position',[1 1 7.8 4.75]);
    ax=axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    ymax=max(P.MeanValidatedAMI+P.StdAcrossReplicates);
    ymin=min(P.MeanValidatedAMI-P.StdAcrossReplicates);
    span=max(ymax-ymin,0.1);
    yLimits=[ymin-0.07*span ymax+0.10*span];

    % M32 is infeasible beyond d_min=2/(M-1); mark the boundary without
    % extending the M32 line into the infeasible region.
    d32max=2*result.config.P_avg/(32-1);
    patch(ax,[d32max 0.104 0.104 d32max], ...
        [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
        [0.78 0.78 0.78], ...
        'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    xline(ax,d32max,':','Color',[0.45 0.45 0.45], ...
        'LineWidth',1,'HandleVisibility','off');

    colors=[0 0.4470 0.7410;0.8500 0.3250 0.0980];
    markers={'o','s'};
    for iM=1:2
        m=[16 32];
        R=P(P.M==m(iM),:);
        errorbar(ax,R.MinGap,R.MeanValidatedAMI,R.StdAcrossReplicates, ...
            'LineStyle','-','Color',colors(iM,:), ...
            'LineWidth',1.5,'Marker',markers{iM}, ...
            'MarkerSize',6,'MarkerFaceColor','w','CapSize',5, ...
            'DisplayName',sprintf('M = %d',m(iM)));
    end

    xlabel(ax,'Minimum spacing d_{min}','Interpreter','tex');
    ylabel(ax,'Validated AMI [bits/symbol]');
    title(ax,'SNR = 20 dB, \sigma_R^2 = 0.2', ...
        'Interpreter','tex','FontWeight','normal');
    legend(ax,'Location','northeast');

    xticks(ax,[0 0.01 0.02 0.03 0.04 0.05 0.07 0.10]);
    xlim(ax,[-0.003 0.104]);
    ylim(ax,yLimits);

    ax.Color='w';
    ax.XColor='k';ax.YColor='k';
    ax.FontName='Times New Roman';ax.FontSize=10;ax.LineWidth=0.8;
    ax.GridColor=[0.82 0.82 0.82];ax.GridAlpha=0.5;
    try,ax.Toolbar.Visible='off';catch,end
    drawnow;

    stem='Fig_K4_combined_min_gap';
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

    fprintf('Combined spacing figure saved:\n  %s\n  %s\n  %s\n', ...
        paths.pdf,paths.png,paths.fig);

    if o.CloseFigure,close(fig);end
end


function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
