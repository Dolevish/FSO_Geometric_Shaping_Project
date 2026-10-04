function out = plot_i13d_optimizer_benchmark(result,varargin)
%PLOT_I13D_OPTIMIZER_BENCHMARK  Reviewer-facing I.13D comparison figures.
%
% Produces:
%   1) Fig_I13D_SA_vs_PatternSearch
%      (a) mean validated AMI: SA vs Pattern Search with y=x
%      (b) delta I = I_SA-I_PS with +/-1e-3 equivalence band.
%
%   2) Fig_I13D_four_method_comparison (when Zhou I.13B data are present)
%      Four panels by M comparing PAM, Zhou-2025 GS, Pattern Search GS and
%      frozen-SA GS at the four benchmark channel conditions.
%
%   3) Fig_I13D_restart_variability
%      Internal diagnostic showing every Pattern Search restart and the
%      frozen SA winner for each optimization case.
%
% The function only plots saved benchmark results; it performs no AMI
% evaluation or optimization.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'CloseFigures',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    if ischar(result)||isstring(result)
        S=load(char(result),'result');
        result=S.result;
    end
    if ~isstruct(result) || ~isfield(result,'physicalSummary')
        error('plot_i13d_optimizer_benchmark:Input','Expected I.13D result struct or MAT path.');
    end

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=result.outputDirectory;
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    P=sortrows(result.physicalSummary,{'M','SNRdB','SigmaR2'});
    C=sortrows(result.caseSummary,{'M','SNRdB','SigmaR2','Replicate'});
    R=sortrows(result.restartSummary,{'M','SNRdB','SigmaR2','Replicate','Restart'});
    tol=result.config.tieToleranceBits;

    f1=make_primary(P,tol);
    out.primary=save_triplet(f1,outDir,'Fig_I13D_SA_vs_PatternSearch');

    if ismember('ZhouValidated',P.Properties.VariableNames) && all(isfinite(P.ZhouValidated))
        f2=make_four_method(P);
        out.fourMethod=save_triplet(f2,outDir,'Fig_I13D_four_method_comparison');
    else
        f2=[];
        out.fourMethod=struct();
    end

    f3=make_restart_variability(R,C);
    out.restartVariability=save_triplet(f3,outDir,'Fig_I13D_restart_variability');

    if o.CloseFigures
        close(f1);
        if ~isempty(f2),close(f2);end
        close(f3);
    end
end


function fig=make_primary(P,tol)
    MVec=unique(P.M).';
    colors=lines(numel(MVec));
    markers={'o','s'};

    fig=figure('Color','w','Name','I13D SA vs Pattern Search', ...
        'Units','inches','Position',[1 1 10.5 4.4]);
    tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    allVals=[P.SAMean;P.PSMean];
    lo=min(allVals); hi=max(allVals);
    pad=max(0.03,0.04*(hi-lo+eps));
    plot(ax,[lo-pad hi+pad],[lo-pad hi+pad],'k--','LineWidth',1.2, ...
        'DisplayName','y=x');

    legendHandles=gobjects(numel(MVec),1);
    for iM=1:numel(MVec)
        M=MVec(iM);
        rows=find(P.M==M);
        for q=1:numel(rows)
            k=rows(q);
            iSig=1+(P.SigmaR2(k)>0);
            face='none';
            if P.SNRdB(k)>=30,face=colors(iM,:);end
            h=scatter(ax,P.SAMean(k),P.PSMean(k),55, ...
                'Marker',markers{iSig},'MarkerEdgeColor',colors(iM,:), ...
                'MarkerFaceColor',face,'LineWidth',1.2);
            if q==1,legendHandles(iM)=h;end
        end
    end
    xlabel(ax,'SA validated AMI [bits/symbol]');
    ylabel(ax,'Pattern Search validated AMI [bits/symbol]');
    title(ax,'(a) Optimizer agreement');
    xlim(ax,[lo-pad hi+pad]); ylim(ax,[lo-pad hi+pad]);
    legLabels=arrayfun(@(m)sprintf('M=%d',m),MVec,'UniformOutput',false);
    legend(ax,legendHandles,legLabels,'Location','best','FontSize',8);
    set(ax,'FontSize',9,'LineWidth',0.8);

    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    n=height(P); x=1:n;
    yl=max([max(abs(P.SAMinusPSMean))*1.15,5*tol,0.005]);

    patch(ax,[0.5 n+0.5 n+0.5 0.5],[-tol -tol tol tol], ...
        [0.85 0.85 0.85],'FaceAlpha',0.35,'EdgeColor','none', ...
        'DisplayName',sprintf('|\\Delta I| \\le %.0e',tol));
    yline(ax,0,'k-','LineWidth',1,'HandleVisibility','off');
    yline(ax,tol,'k:','LineWidth',0.8,'HandleVisibility','off');
    yline(ax,-tol,'k:','LineWidth',0.8,'HandleVisibility','off');

    for iM=1:numel(MVec)
        idx=P.M==MVec(iM);
        scatter(ax,x(idx),P.SAMinusPSMean(idx),48, ...
            'MarkerEdgeColor',colors(iM,:),'MarkerFaceColor',colors(iM,:), ...
            'DisplayName',sprintf('M=%d',MVec(iM)));
    end

    labels=cell(n,1);
    for k=1:n
        labels{k}=sprintf('%d|%d|%.1f',P.M(k),P.SNRdB(k),P.SigmaR2(k));
    end
    xlim(ax,[0.5 n+0.5]); ylim(ax,[-yl yl]);
    xticks(ax,x); xticklabels(ax,labels); xtickangle(ax,55);
    ylabel(ax,'I_{SA}-I_{PS} [bits/symbol]');
    xlabel(ax,'M | SNR [dB] | \sigma_R^2');
    title(ax,'(b) Difference with equivalence band');
    legend(ax,'Location','best','FontSize',8);
    set(ax,'FontSize',8.5,'LineWidth',0.8);
end


function fig=make_four_method(P)
    MVec=unique(P.M).';
    fig=figure('Color','w','Name','I13D Four-method comparison', ...
        'Units','inches','Position',[1 1 8.0 6.0]);
    tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    for iM=1:numel(MVec)
        M=MVec(iM);
        R=sortrows(P(P.M==M,:),{'SNRdB','SigmaR2'});
        ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        x=1:height(R);
        plot(ax,x,R.PAMValidated,'--o','LineWidth',1.2,'MarkerSize',4, ...
            'DisplayName','Uniform PAM');
        plot(ax,x,R.ZhouValidated,'-d','LineWidth',1.3,'MarkerSize',4, ...
            'DisplayName','Zhou GS');
        plot(ax,x,R.PSMean,'-s','LineWidth',1.5,'MarkerSize',4, ...
            'DisplayName','Pattern Search');
        plot(ax,x,R.SAMean,'-^','LineWidth',1.6,'MarkerSize',4, ...
            'DisplayName','SA');

        labels=cell(height(R),1);
        for k=1:height(R)
            labels{k}=sprintf('%d/%.1f',R.SNRdB(k),R.SigmaR2(k));
        end
        xticks(ax,x); xticklabels(ax,labels);
        xlabel(ax,'SNR / \sigma_R^2');
        ylabel(ax,'Validated AMI [bits/symbol]');
        title(ax,sprintf('M = %d',M),'FontWeight','normal');
        if iM==1,legend(ax,'Location','best','FontSize',7);end
        set(ax,'FontSize',8.5,'LineWidth',0.8);
    end
end


function fig=make_restart_variability(R,C)
    MVec=unique(C.M).';
    fig=figure('Color','w','Name','I13D restart variability', ...
        'Units','inches','Position',[1 1 9.0 6.2]);
    tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    for iM=1:numel(MVec)
        M=MVec(iM);
        cases=sortrows(C(C.M==M,:),{'SNRdB','SigmaR2','Replicate'});
        ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        for c=1:height(cases)
            Q=R(R.M==M & R.SNRdB==cases.SNRdB(c) & ...
                abs(R.SigmaR2-cases.SigmaR2(c))<1e-12 & ...
                R.Replicate==cases.Replicate(c),:);
            x=c+linspace(-0.16,0.16,height(Q)).';
            scatter(ax,x,Q.PSValidated,16,'filled','MarkerFaceAlpha',0.65, ...
                'HandleVisibility','off');
            plot(ax,c,cases.SAValidated(c),'kx','MarkerSize',7,'LineWidth',1.3, ...
                'HandleVisibility','off');
        end

        labels=cell(height(cases),1);
        for c=1:height(cases)
            labels{c}=sprintf('%d/%.1f/r%d',cases.SNRdB(c),cases.SigmaR2(c),cases.Replicate(c));
        end
        xticks(ax,1:height(cases)); xticklabels(ax,labels); xtickangle(ax,55);
        ylabel(ax,'Validated AMI [bits/symbol]');
        xlabel(ax,'SNR / \sigma_R^2 / replicate');
        title(ax,sprintf('M = %d: PS restarts (dots), SA winner (x)',M), ...
            'FontWeight','normal');
        set(ax,'FontSize',8,'LineWidth',0.8);
    end
end


function p=save_triplet(fig,outDir,stem)
    p=struct();
    p.png=fullfile(outDir,[stem '.png']);
    p.pdf=fullfile(outDir,[stem '.pdf']);
    p.fig=fullfile(outDir,[stem '.fig']);

    if ~isgraphics(fig,'figure')
        error('plot_i13d_optimizer_benchmark:InvalidFigure', ...
            'Expected a valid figure handle before export.');
    end

    ax=findall(fig,'Type','axes');
    for k=1:numel(ax)
        try
            ax(k).Toolbar.Visible='off';
        catch
        end
    end
    drawnow;

    try
        savefig(fig,p.fig);
    catch ME1
        if isgraphics(fig,'figure')
            try
                hgsave(fig,p.fig);
            catch ME2
                error('plot_i13d_optimizer_benchmark:FigSaveFailed', ...
                    'Could not save FIG (%s; fallback: %s).',ME1.message,ME2.message);
            end
        else
            error('plot_i13d_optimizer_benchmark:FigureInvalidated', ...
                'Figure handle became invalid before FIG save: %s',ME1.message);
        end
    end

    exportgraphics(fig,p.png,'Resolution',300);
    try
        exportgraphics(fig,p.pdf,'ContentType','vector');
    catch
        exportgraphics(fig,p.pdf);
    end
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
