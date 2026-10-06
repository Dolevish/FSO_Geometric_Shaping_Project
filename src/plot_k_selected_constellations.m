function out = plot_k_selected_constellations(resultOrPath,varargin)
%PLOT_K_SELECTED_CONSTELLATIONS  Commit-K four-panel constellation figures.
%
% Creates one publication-style figure for each constellation order M.
% Each figure contains four panels for
%       d_min = [0.005 0.01 0.02 0.05].
%
% In every panel, the selected d_min-constrained constellation is compared
% against the selected d_min=0 constellation for the same M. Selection is
% performed independently at every (M,d_min): choose the replicate with the
% highest validated AMI (deterministic tie-break: lower replicate number).
%
% Visual style follows the project's legacy constellation-comparison figure:
% intensity on the x-axis, two stem/lollipop series at different display
% heights, a legend, a central condition box, and AMI annotations.
%
% No SA or AMI evaluation is rerun.
%
% Usage:
%   Fk = plot_k_selected_constellations(Rk);
%
% or after restarting MATLAB:
%   Fk = plot_k_selected_constellations('/path/to/K_min_gap_results.mat');
%
% Optional name-value arguments:
%   'OutputDirectory'  default <K run>/selected_constellation_figures
%   'MinGapsToPlot'    default [0.005 0.01 0.02 0.05]
%   'CloseFigures'     default false
%
% Outputs per M:
%   Fig_K_selected_constellations_M16.{pdf,png,fig}
%   Fig_K_selected_constellations_M32.{pdf,png,fig}
% plus K_selected_constellation_plot_summary.csv.
%
% No global-optimum claim is implied.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'MinGapsToPlot',[0.005 0.01 0.02 0.05],@valid_gap_vector);
    addParameter(p,'CloseFigures',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    result=load_result(resultOrPath);
    validate_result(result);

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    elseif isfield(result,'outputDirectory') && ~isempty(result.outputDirectory)
        outDir=fullfile(result.outputDirectory,'selected_constellation_figures');
    else
        outDir=fullfile(fso_result_utils.revision2_results_root(), ...
            'commit_k_min_gap','selected_constellation_figures');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('plot_k_selected_constellations:mkdir','%s',msg);end
    end

    selected=select_best_replicates(result.replicateSummary);
    dPlot=double(o.MinGapsToPlot(:).');
    MVec=unique(selected.M).';

    if ~isequal(sort(dPlot),sort([0.005 0.01 0.02 0.05]))
        warning('plot_k_selected_constellations:NonCanonicalGapSet', ...
            'The canonical Commit-K four-panel set is [0.005 0.01 0.02 0.05].');
    end

    paths=struct();
    figs=struct();

    summaryRows=cell(numel(MVec)*numel(dPlot),1);
    q=0;

    fprintf('\n============================================================\n');
    fprintf('COMMIT K - SELECTED CONSTELLATION FIGURES\n');
    fprintf('Selection: highest validated AMI replicate at each (M,d_min).\n');
    fprintf('Reference: selected d_min=0 constellation for the same M.\n');
    fprintf('Panels: d_min=%s\n',mat2str(dPlot));
    fprintf('Output: %s\n',outDir);
    fprintf('No optimization or AMI evaluation is rerun.\n');
    fprintf('============================================================\n');

    for iM=1:numel(MVec)
        M=MVec(iM);
        rowsM=selected(selected.M==M,:);

        ref=rowsM(abs(rowsM.MinGap)<1e-12,:);
        if height(ref)~=1
            error('plot_k_selected_constellations:Reference', ...
                'Expected exactly one selected d_min=0 row for M=%d.',M);
        end
        xRef=sort(double(ref.WinnerConstellation{1}(:)));
        refAMI=double(ref.WinnerValidatedAMI);
        refRep=double(ref.Replicate);

        % Use one common x-axis across all four panels. This makes the
        % geometry changes directly comparable and prevents per-panel
        % autoscaling from exaggerating or hiding movement of high levels.
        allX=xRef(:);
        for d=dPlot
            row=rowsM(abs(rowsM.MinGap-d)<1e-12,:);
            if height(row)~=1
                error('plot_k_selected_constellations:GapRow', ...
                    'Expected one selected row for M=%d, d_min=%.6g.',M,d);
            end
            allX=[allX; double(row.WinnerConstellation{1}(:))]; %#ok<AGROW>
        end
        [xLim,xTicks]=shared_intensity_axis(allX);

        fig=figure('Name',sprintf('Commit K selected constellations M=%d',M), ...
            'Color','w','Units','inches','Position',[0.4 0.4 11.4 6.7]);
        tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('M = %d, SNR = 20 dB, \\sigma_R^2 = 0.2',M), ...
            'Interpreter','tex','FontName','Times New Roman', ...
            'FontSize',15,'FontWeight','normal');

        for j=1:numel(dPlot)
            d=dPlot(j);
            row=rowsM(abs(rowsM.MinGap-d)<1e-12,:);
            xSel=sort(double(row.WinnerConstellation{1}(:)));
            selAMI=double(row.WinnerValidatedAMI);
            selRep=double(row.Replicate);
            dActual=min(diff(xSel));
            deltaAMI=selAMI-refAMI;

            ax=nexttile(tl,j);
            hold(ax,'on'); grid(ax,'on'); box(ax,'on');

            hRef=stem(ax,xRef,ones(size(xRef)),'o', ...
                'Color',[0 0.4470 0.7410], ...
                'MarkerSize',marker_size(M),'LineWidth',1.25, ...
                'MarkerFaceColor','w','DisplayName','GS reference, d_{min}=0');

            hSel=stem(ax,xSel,0.52*ones(size(xSel)),'x', ...
                'Color',[0.8500 0.3250 0.0980], ...
                'MarkerSize',marker_size(M)+1,'LineWidth',1.55, ...
                'DisplayName',sprintf('GS, d_{min}=%.3g',d));

            xlim(ax,xLim);
            xticks(ax,xTicks);
            ylim(ax,[-0.08 1.34]);
            yticks(ax,[0.52 1.0]);
            yticklabels(ax,{'',''});
            xlabel(ax,'Intensity level','FontName','Times New Roman', ...
                'FontSize',10,'FontWeight','bold');

            % Condition box in the same spirit as the supplied legacy figure.
            text(ax,0.50,0.965,sprintf(['d_{min}=%.3g\n' ...
                'd_{actual}=%.4f\n\\Delta AMI=%+.4f bits'], ...
                d,dActual,deltaAMI), ...
                'Units','normalized','HorizontalAlignment','center', ...
                'VerticalAlignment','top','FontName','Times New Roman', ...
                'FontSize',8.5,'FontWeight','bold','Interpreter','tex', ...
                'BackgroundColor','white','EdgeColor','black','Margin',4);

            % AMI annotations at top right.
            text(ax,0.975,0.97,sprintf('REF: %.4f',refAMI), ...
                'Units','normalized','HorizontalAlignment','right', ...
                'VerticalAlignment','top','FontName','Times New Roman', ...
                'FontSize',8.5,'FontWeight','bold','Color',[0 0.4470 0.7410]);
            text(ax,0.975,0.855,sprintf('OPT: %.4f',selAMI), ...
                'Units','normalized','HorizontalAlignment','right', ...
                'VerticalAlignment','top','FontName','Times New Roman', ...
                'FontSize',8.5,'FontWeight','bold','Color',[0.8500 0.3250 0.0980]);

            legend(ax,[hRef hSel], ...
                {'GS reference, d_{min}=0',sprintf('GS, d_{min}=%.3g',d)}, ...
                'Location','northwest','FontName','Times New Roman', ...
                'FontSize',7.5,'Interpreter','tex','Box','on');

            ax.FontName='Times New Roman';
            ax.FontSize=8.5;
            ax.LineWidth=0.8;
            try,ax.Toolbar.Visible='off';catch,end

            q=q+1;
            summaryRows{q}=struct( ...
                'M',M,'MinGap',d,'SelectedReplicate',selRep, ...
                'ReferenceReplicate',refRep,'ValidatedAMI',selAMI, ...
                'ReferenceValidatedAMI',refAMI,'DeltaAMI',deltaAMI, ...
                'ActualMinGap',dActual,'MaxLevel',max(xSel));
        end

        stemName=sprintf('Fig_K_selected_constellations_M%d',M);
        P=save_triplet(fig,outDir,stemName);
        key=sprintf('M%d',M);
        paths.(key)=P;
        figs.(key)=fig;

        fprintf('M=%d figure saved:\n  %s\n  %s\n  %s\n', ...
            M,P.pdf,P.png,P.fig);
    end

    summaryRows=summaryRows(1:q);
    plotSummary=struct2table(vertcat(summaryRows{:}));
    plotSummary=sortrows(plotSummary,{'M','MinGap'});
    writetable(plotSummary,fullfile(outDir,'K_selected_constellation_plot_summary.csv'));

    out=struct();
    out.outputDirectory=outDir;
    out.paths=paths;
    out.figures=figs;
    out.selectedTable=selected;
    out.plotSummary=plotSummary;
    out.referencePolicy='best validated replicate at d_min=0 for each M';
    out.selectionPolicy='best validated replicate independently at each (M,d_min)';
    out.globalOptimumClaim=false;

    fprintf('Summary CSV: %s\n', ...
        fullfile(outDir,'K_selected_constellation_plot_summary.csv'));
    fprintf('============================================================\n');

    if o.CloseFigures
        names=fieldnames(figs);
        for k=1:numel(names)
            if isgraphics(figs.(names{k}))
                close(figs.(names{k}));
            end
        end
    end
end


function selected=select_best_replicates(R)
    keys=unique(R(:,{'M','SNRdB','SigmaR2','MinGap'}),'rows','sorted');
    pieces=cell(height(keys),1);

    for k=1:height(keys)
        idx=R.M==keys.M(k) & ...
            abs(R.SNRdB-keys.SNRdB(k))<1e-12 & ...
            abs(R.SigmaR2-keys.SigmaR2(k))<1e-12 & ...
            abs(R.MinGap-keys.MinGap(k))<1e-12;
        G=R(idx,:);
        if isempty(G)
            error('plot_k_selected_constellations:MissingGroup', ...
                'Missing replicate group at row %d.',k);
        end
        G=sortrows(G,{'WinnerValidatedAMI','Replicate'},{'descend','ascend'});
        pieces{k}=G(1,:);
    end

    selected=vertcat(pieces{:});
    selected=sortrows(selected,{'M','MinGap'});
end


function [lim,ticks]=shared_intensity_axis(x)
    x=double(x(:));
    x=x(isfinite(x));
    if isempty(x)
        error('plot_k_selected_constellations:EmptyGeometry','No finite intensity levels.');
    end

    xmin=min(x); xmax=max(x);
    span=max(xmax-xmin,1);
    margin=max(0.04*span,0.08);

    lo=min(-margin,xmin-margin);
    rawHi=xmax+margin;

    % Nice common tick spacing based on full intensity span.
    targetTicks=7;
    rawStep=max((rawHi-lo)/targetTicks,eps);
    decade=10^floor(log10(rawStep));
    scaled=rawStep/decade;
    if scaled<=1
        nice=1;
    elseif scaled<=2
        nice=2;
    elseif scaled<=2.5
        nice=2.5;
    elseif scaled<=5
        nice=5;
    else
        nice=10;
    end
    step=nice*decade;

    tickStart=floor(max(0,lo)/step)*step;
    tickEnd=ceil(rawHi/step)*step;
    ticks=tickStart:step:tickEnd;
    if isempty(ticks),ticks=[0 tickEnd];end
    lim=[lo tickEnd+0.02*max(step,span)];
end


function s=marker_size(M)
    if M>=32,s=4.5;
    elseif M>=16,s=5.5;
    else,s=7;
    end
end


function P=save_triplet(fig,outDir,stem)
    P=struct();
    P.pdf=fullfile(outDir,[stem '.pdf']);
    P.png=fullfile(outDir,[stem '.png']);
    P.fig=fullfile(outDir,[stem '.fig']);

    try,savefig(fig,P.fig);catch,hgsave(fig,P.fig);end
    exportgraphics(fig,P.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,P.pdf,'ContentType','vector','BackgroundColor','white');
    catch
        exportgraphics(fig,P.pdf,'BackgroundColor','white');
    end
end


function validate_result(result)
    if ~isstruct(result)||~isfield(result,'replicateSummary')|| ...
            ~istable(result.replicateSummary)
        error('plot_k_selected_constellations:Input', ...
            'Expected Commit-K result struct or MAT path with replicateSummary.');
    end

    R=result.replicateSummary;
    required={'M','SNRdB','SigmaR2','MinGap','Replicate', ...
        'WinnerValidatedAMI','WinnerConstellation'};
    missing=setdiff(required,R.Properties.VariableNames);
    if ~isempty(missing)
        error('plot_k_selected_constellations:Fields', ...
            'replicateSummary is missing: %s',strjoin(missing,', '));
    end

    if ~all(ismember([16 32],unique(R.M).'))
        error('plot_k_selected_constellations:MGrid', ...
            'Commit-K M=16 and M=32 rows are required.');
    end
    expected=[0 0.005 0.01 0.02 0.05];
    actual=unique(R.MinGap).';
    if numel(actual)<numel(expected) || any(~ismember(expected,actual))
        error('plot_k_selected_constellations:GapGrid', ...
            'Commit-K d_min grid is incomplete.');
    end
end


function result=load_result(v)
    if isstruct(v),result=v;return;end
    if ~(ischar(v)||(isstring(v)&&isscalar(v)))
        error('plot_k_selected_constellations:InputType', ...
            'Input must be a Commit-K result struct or MAT-file path.');
    end
    path=char(string(v));
    if exist(path,'file')~=2
        error('plot_k_selected_constellations:MissingFile','Missing MAT file: %s',path);
    end
    S=load(path);
    if isfield(S,'result'),result=S.result;
    elseif isfield(S,'Rk'),result=S.Rk;
    else,error('plot_k_selected_constellations:Variable', ...
            'MAT file must contain variable result or Rk.');
    end
end


function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=valid_gap_vector(v)
    tf=isnumeric(v)&&isvector(v)&&numel(v)==4&&all(isfinite(v))&&all(v>0);
end
