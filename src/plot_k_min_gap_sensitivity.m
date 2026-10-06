function paths = plot_k_min_gap_sensitivity(result,varargin)
%PLOT_K_MIN_GAP_SENSITIVITY  One-axis Commit-K AMI-vs-d_min figure.
%
% Each marker is the mean of two independently validated replicate winners.
% Error bars are the sample standard deviation across the two replicates.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'CloseFigure',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    if ischar(result)||isstring(result)
        S=load(char(result));
        if isfield(S,'result'),result=S.result;
        else,error('plot_k_min_gap_sensitivity:MAT','MAT file must contain result.');end
    end
    if ~isstruct(result)||~isfield(result,'physicalSummary')
        error('plot_k_min_gap_sensitivity:Input','Expected Commit-K result struct or MAT path.');
    end

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=result.outputDirectory;
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    T=sortrows(result.physicalSummary,{'M','MinGap'});
    MVec=unique(T.M).';

    fig=figure('Color','w','Name','Commit K minimum-spacing sensitivity', ...
        'Units','inches','Position',[1 1 6.9 4.25]);
    ax=axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    markers={'o','s'};
    for iM=1:numel(MVec)
        R=T(T.M==MVec(iM),:);
        errorbar(ax,R.MinGap,R.MeanValidatedAMI,R.StdAcrossReplicates, ...
            'LineWidth',1.35,'Marker',markers{iM},'MarkerSize',6, ...
            'MarkerFaceColor','w','CapSize',5, ...
            'DisplayName',sprintf('M = %d',MVec(iM)));
    end

    xlabel(ax,'Minimum spacing d_{min}','Interpreter','tex');
    ylabel(ax,'Validated AMI [bits/symbol]');
    title(ax,'SNR = 20 dB, \sigma_R^2 = 0.2','FontWeight','normal','Interpreter','tex');
    legend(ax,'Location','best');
    xticks(ax,result.config.MinGapVec);
    xlim(ax,[min(result.config.MinGapVec)-0.001 max(result.config.MinGapVec)+0.003]);

    ax.Color='w';
    ax.XColor='k'; ax.YColor='k';
    ax.FontName='Times New Roman'; ax.FontSize=9; ax.LineWidth=0.8;
    ax.GridColor=[0.82 0.82 0.82]; ax.GridAlpha=0.55;
    try,ax.Toolbar.Visible='off';catch,end
    drawnow;

    stem='Fig_K_min_gap_sensitivity';
    paths=struct();
    paths.png=fullfile(outDir,[stem '.png']);
    paths.pdf=fullfile(outDir,[stem '.pdf']);
    paths.fig=fullfile(outDir,[stem '.fig']);

    try,savefig(fig,paths.fig);catch,hgsave(fig,paths.fig);end
    exportgraphics(fig,paths.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,paths.pdf,'ContentType','vector','BackgroundColor','white');
    catch
        exportgraphics(fig,paths.pdf,'BackgroundColor','white');
    end

    if o.CloseFigure,close(fig);end
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
