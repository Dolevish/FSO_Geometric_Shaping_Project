function out = make_ieee_revision_figures_legacy_style(varargin)
%MAKE_IEEE_REVISION_FIGURES_LEGACY_STYLE
% I.12C: regenerate revised IEEE figures using the visual format of the
% pre-revision manuscript figures, while using only the frozen I.11.1
% production results.
%
% The old manuscript figure layout is reproduced:
%   1) AMI vs SNR for M=8:
%      turbulence encoded by color/marker, GS solid, PAM dashed.
%   2) M=8, SNR=20 dB constellation comparison:
%      four turbulence panels, intensity on x-axis, PAM/GS stem levels.
%   3) AMI vs SNR for all M:
%      four turbulence panels, M encoded by color/marker, GS solid,
%      PAM dashed.
%
% IMPORTANT SCIENTIFIC POLICIES:
%   - No SA or AMI evaluation is performed here.
%   - Performance curves use the mean validated AMI of the two production
%     replicates, as frozen in I.12A.
%   - The constellation geometry uses replicate 1 consistently, as frozen
%     in I.12A-v2. The AMI annotation beside that geometry therefore uses
%     replicate 1, not the two-replicate mean.
%   - The implementation parameter is sigma_X^2 (normalized intensity
%     variance/scintillation index), not Rytov variance. By default the
%     corrected symbol sigma_X^2 is shown. Set 'UseLegacySigmaSymbol',true
%     only when an exact old-label rendering is desired.
%   - The legacy axis policy is retained. Where new production values exceed
%     an old panel limit, the upper limit is expanded just enough to avoid
%     clipping while preserving the old tick spacing.
%
% Usage:
%   out = make_ieee_revision_figures_legacy_style();
%   out = make_ieee_revision_figures_legacy_style( ...
%       'ProductionBundlePath','/.../i11p1ProductionBundle.mat');
%
% Outputs are written under:
%   <production>/i12_analysis/ieee_figures_legacy_style/
%
% No global-optimum claim is implied.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'UseLegacySigmaSymbol',false,@logical_scalar);
    addParameter(p,'CloseFigures',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    A=analyze_i11p1_production_results( ...
        'ProductionBundlePath',o.ProductionBundlePath, ...
        'ResultsRoot',o.ResultsRoot,'WriteCSV',true);

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=fullfile(A.outputDirectory,'ieee_figures_legacy_style');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('make_ieee_revision_figures_legacy_style:mkdir','%s',msg);end
    end

    if o.UseLegacySigmaSymbol
        sigmaSymbol='\\sigma_R^2';
    else
        sigmaSymbol='\\sigma_X^2';
    end

    P=A.physicalTable;
    C=A.fig4_M8_SNR20;
    validate_inputs(P,C);

    fAMI=make_m8_ami_legacy(P,sigmaSymbol);
    fConst=make_constellation_legacy(C,A.tableI_M8_SNR20,sigmaSymbol);
    fAllM=make_allm_legacy(P,sigmaSymbol);

    out=struct();
    out.version='I12C-v1';
    out.sourceBundlePath=A.sourceBundlePath;
    out.outputDirectory=outDir;
    out.sigmaSymbol=sigmaSymbol;
    out.axisPolicy='legacy-paper axes; expand upper bound only to prevent clipping';
    out.performanceAggregation='mean validated AMI across two replicates';
    out.geometryPolicy='replicate 1';
    out.amiVsSNR=save_triplet(fAMI,outDir,'AMI vs SNR');
    out.constellation=save_triplet(fConst,outDir,'Constellation Comparison2');
    out.amiVsSNRAllM=save_triplet(fAllM,outDir,'AMI vs SNR allM');

    save(fullfile(outDir,'i12c_legacy_figure_manifest.mat'),'out','-v7.3');

    fprintf('\n============================================================\n');
    fprintf('I.12C LEGACY-STYLE IEEE FIGURES COMPLETE\n');
    fprintf('Source bundle : %s\n',A.sourceBundlePath);
    fprintf('Sigma symbol  : %s\n',sigmaSymbol);
    fprintf('Output dir    : %s\n',outDir);
    fprintf('AMI vs SNR    : %s\n',out.amiVsSNR.pdf);
    fprintf('Constellation : %s\n',out.constellation.pdf);
    fprintf('AMI all M     : %s\n',out.amiVsSNRAllM.pdf);
    fprintf('============================================================\n');

    if o.CloseFigures
        close(fAMI); close(fConst); close(fAllM);
    end
end


function fig=make_m8_ami_legacy(P,sigmaSymbol)
% Match the original manuscript's M=8 figure styling.

    M=8;
    sigVec=[0 0.1 0.2 0.3];
    snrVec=[5 10 15 20 25 30];

    colors={ ...
        [0.00 0.45 0.74], ... % blue
        [0.64 0.08 0.18], ... % dark red
        [0.47 0.67 0.19], ... % green
        [0.93 0.69 0.13]  ... % gold
    };
    markers={'s','*','o','d'};

    fig=figure('Name','AMI vs SNR','Color','w', ...
        'Units','normalized','Position',[0.02 0.52 0.45 0.40]);
    ax=axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    for i=1:numel(sigVec)
        R=sortrows(P(P.M==M & abs(P.SigmaX2-sigVec(i))<1e-12,:),'SNRdB');
        assert(height(R)==numel(snrVec),'Incomplete M=8 AMI curve.');
        c=colors{i}; m=markers{i};

        plot(ax,R.SNRdB,R.GSMeanValidated,['-' m], ...
            'Color',c,'LineWidth',2.0, ...
            'MarkerFaceColor',c,'MarkerSize',7, ...
            'DisplayName',sprintf('GS opt  (%s=%.2f)',sigmaSymbol,sigVec(i)));

        plot(ax,R.SNRdB,R.PAMValidated,['--' m], ...
            'Color',c,'LineWidth',1.5, ...
            'MarkerFaceColor','w','MarkerSize',6, ...
            'DisplayName',sprintf('PAM     (%s=%.2f)',sigmaSymbol,sigVec(i)));
    end

    yline(ax,log2(M),'k:','LineWidth',1,'HandleVisibility','off');
    text(ax,snrVec(end)+0.30,log2(M),sprintf('log_2(%d)=%d',M,log2(M)), ...
        'FontSize',9,'VerticalAlignment','bottom');

    xlabel(ax,'SNR [dB]','FontSize',12,'FontWeight','bold');
    ylabel(ax,'AMI [bits/symbol]','FontSize',12,'FontWeight','bold');
    legend(ax,'Location','southeast','NumColumns',2,'FontSize',10);

    xlim(ax,[snrVec(1)-1 snrVec(end)+1]);
    xticks(ax,snrVec);
    ylim(ax,[0 log2(M)*1.08]);
    yticks(ax,0:0.5:3.0);
    set(ax,'FontSize',11,'LineWidth',0.8);
end


function fig=make_constellation_legacy(C,tableI,sigmaSymbol)
% Match the original manuscript's four-panel stem constellation figure.

    C=sortrows(C,'SigmaX2');
    tableI=sortrows(tableI,'SigmaX2');
    sigVec=C.SigmaX2(:).';
    M=8;

    if ~isequal(sigVec,[0 0.1 0.2 0.3])
        error('make_ieee_revision_figures_legacy_style:ConstellationGrid', ...
            'Expected sigma grid [0 0.1 0.2 0.3].');
    end

    % Original-paper target axis ranges/ticks. The right edge is expanded
    % when the revised constellation would otherwise be clipped.
    oldXMax=[2.5 4.0 3.0 3.5];
    xTickStep=[0.5 1.0 0.5 0.5];

    fig=figure('Name','Constellation Comparison','Color','w', ...
        'Units','normalized','Position',[0.02 0.15 0.96 0.75]);

    for i=1:numel(sigVec)
        ax=subplot(2,2,i,'Parent',fig);
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        xPam=zeros(1,M); xGS=zeros(1,M);
        for j=1:M
            xPam(j)=C.(sprintf('PAM_x%d',j))(i);
            xGS(j)=C.(sprintf('GS_representative_x%d',j))(i);
        end
        xPam=sort(xPam); xGS=sort(xGS);

        hPam=stem(ax,xPam,ones(size(xPam)),'bo', ...
            'MarkerSize',8,'LineWidth',1.5,'MarkerFaceColor','none', ...
            'DisplayName','Uniform PAM');
        hGS=stem(ax,xGS,0.5*ones(size(xGS)),'rx', ...
            'MarkerSize',10,'LineWidth',2.0, ...
            'DisplayName','GS Optimized');

        ylim(ax,[-0.1 1.4]);
        yticks(ax,[0.5 1.0]);
        yticklabels(ax,{'',''});

        step=xTickStep(i);
        targetMax=oldXMax(i);
        dataMax=max([xPam xGS]);
        xMax=max(targetMax,ceil((dataMax+0.05)/step)*step);
        xMin=-0.2;
        if i==2,xMin=-0.4;end
        xlim(ax,[xMin xMax+0.02]);
        xticks(ax,0:step:xMax);

        xlabel(ax,'Intensity level','FontSize',12,'FontWeight','bold');

        pamAMI=tableI.PAMValidated(i);
        gsRep1AMI=C.GSRep1AMI(i);
        gain=gsRep1AMI-pamAMI;

        % Central box matches the paper: sigma on first line, gain below.
        text(ax,0.50,0.96,sprintf('%s = %.1f\nGain = %+0.3f bits', ...
            sigmaSymbol,sigVec(i),gain), ...
            'Units','normalized','HorizontalAlignment','center', ...
            'VerticalAlignment','top','FontSize',10,'FontWeight','bold', ...
            'BackgroundColor','white','EdgeColor','black','Margin',5);

        % Top-right AMI annotations.
        text(ax,0.95,0.97,sprintf('PAM: %.3f',pamAMI), ...
            'Units','normalized','HorizontalAlignment','right', ...
            'VerticalAlignment','top','FontSize',10,'FontWeight','bold', ...
            'Color','b');
        text(ax,0.95,0.85,sprintf('OPT: %.3f',gsRep1AMI), ...
            'Units','normalized','HorizontalAlignment','right', ...
            'VerticalAlignment','top','FontSize',10,'FontWeight','bold', ...
            'Color','r');

        legend(ax,[hPam hGS],{'Uniform PAM','GS Optimized'}, ...
            'Location','northwest','FontSize',9,'Box','on');
        set(ax,'FontSize',10,'LineWidth',0.8);
    end
end


function fig=make_allm_legacy(P,sigmaSymbol)
% Match the old manuscript's all-M figure: one panel per turbulence level.

    MVec=[4 8 16 32];
    sigVec=[0 0.1 0.2 0.3];
    snrVec=[5 10 15 20 25 30];

    colorsM=[ ...
        0.00 0.45 0.74; ... % M=4 blue
        0.85 0.33 0.10; ... % M=8 red-orange
        0.47 0.67 0.19; ... % M=16 green
        0.49 0.18 0.56  ... % M=32 purple
    ];
    markers={'o','s','d','^'};
    linewidthOpt=2.0;
    linewidthPam=1.5;
    markerSize=7;

    % Target limits/ticks read from the pre-revision paper figure.
    % Upper limits expand if revised production data require it.
    oldYMin=[0 0 0 0.25];
    oldYMax=[4.2 2.0 1.65 1.50];
    yStep=[1.0 0.5 0.5 0.5];

    fig=figure('Name','AMI vs SNR allM','Color','w', ...
        'Position',[100 100 1400 900]);

    for iSig=1:numel(sigVec)
        ax=subplot(2,2,iSig,'Parent',fig);
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        hPam=gobjects(1,numel(MVec));
        hOpt=gobjects(1,numel(MVec));
        panelMax=-inf;

        for iM=1:numel(MVec)
            M=MVec(iM);
            R=sortrows(P(P.M==M & abs(P.SigmaX2-sigVec(iSig))<1e-12,:),'SNRdB');
            assert(height(R)==numel(snrVec),'Incomplete all-M AMI curve.');

            hPam(iM)=plot(ax,R.SNRdB,R.PAMValidated,'--', ...
                'Color',colorsM(iM,:),'LineWidth',linewidthPam, ...
                'Marker',markers{iM},'MarkerSize',markerSize-1, ...
                'MarkerFaceColor','none');

            hOpt(iM)=plot(ax,R.SNRdB,R.GSMeanValidated,'-', ...
                'Color',colorsM(iM,:),'LineWidth',linewidthOpt, ...
                'Marker',markers{iM},'MarkerSize',markerSize, ...
                'MarkerFaceColor',colorsM(iM,:));

            panelMax=max(panelMax,max([R.PAMValidated;R.GSMeanValidated]));
        end

        xlabel(ax,'SNR [dB]','FontSize',12,'FontWeight','bold');
        ylabel(ax,'AMI [bits/symbol]','FontSize',12,'FontWeight','bold');
        title(ax,sprintf('%s = %.1f',sigmaSymbol,sigVec(iSig)), ...
            'FontSize',14,'FontWeight','bold');
        xlim(ax,[4 31]);
        xticks(ax,snrVec);

        step=yStep(iSig);
        yMax=max(oldYMax(iSig),ceil((panelMax+0.03)/step)*step);
        ylim(ax,[oldYMin(iSig) yMax]);
        ticks=ceil(oldYMin(iSig)/step)*step:step:yMax;
        yticks(ax,ticks);

        legendEntries=cell(1,2*numel(MVec));
        hAll=gobjects(1,2*numel(MVec));
        for iM=1:numel(MVec)
            legendEntries{2*iM-1}=sprintf('M=%d PAM',MVec(iM));
            legendEntries{2*iM}=sprintf('M=%d OPT',MVec(iM));
            hAll(2*iM-1)=hPam(iM);
            hAll(2*iM)=hOpt(iM);
        end
        legend(ax,hAll,legendEntries,'Location','southeast', ...
            'FontSize',8,'NumColumns',2,'Box','on');

        set(ax,'FontSize',10,'LineWidth',0.8);
    end
end


function paths=save_triplet(fig,outDir,stem)
    paths=struct();
    paths.pdf=fullfile(outDir,[stem '.pdf']);
    paths.png=fullfile(outDir,[stem '.png']);
    paths.fig=fullfile(outDir,[stem '.fig']);

    exportgraphics(fig,paths.png,'Resolution',300);
    try
        exportgraphics(fig,paths.pdf,'ContentType','vector');
    catch
        exportgraphics(fig,paths.pdf);
    end
    savefig(fig,paths.fig);
end


function validate_inputs(P,C)
    if height(P)~=96
        error('make_ieee_revision_figures_legacy_style:PhysicalCount', ...
            'Expected 96 physical points, found %d.',height(P));
    end
    if ~isequal(unique(P.M).',[4 8 16 32]) || ...
            ~isequal(unique(P.SNRdB).',[5 10 15 20 25 30]) || ...
            max(abs(unique(P.SigmaX2).'-[0 0.1 0.2 0.3]))>1e-12
        error('make_ieee_revision_figures_legacy_style:Grid', ...
            'Unexpected production grid.');
    end

    required={'SigmaX2','GSRep1AMI','PAM_x1','GS_representative_x1'};
    missing=setdiff(required,C.Properties.VariableNames);
    if ~isempty(missing)
        error('make_ieee_revision_figures_legacy_style:ConstellationColumns', ...
            'Missing constellation columns: %s',strjoin(missing,', '));
    end
    if height(C)~=4
        error('make_ieee_revision_figures_legacy_style:ConstellationRows', ...
            'Expected four M=8,SNR=20 constellation rows.');
    end
end


function tf=text_scalar(v)
    tf=ischar(v)||(isstring(v)&&isscalar(v));
end

function tf=logical_scalar(v)
    tf=islogical(v)&&isscalar(v);
end
