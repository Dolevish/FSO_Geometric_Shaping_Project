function out = make_i13ef_reviewer12_publication_outputs(result,varargin)
%MAKE_I13EF_REVIEWER12_PUBLICATION_OUTPUTS
% Publication-ready Reviewer-1/2 figures and summary tables from I.13D.
%
% This function performs NO optimization and NO AMI evaluation. It only
% post-processes frozen I.13B/I.13D validated results.
%
% Usage:
%   O = make_i13ef_reviewer12_publication_outputs(Rfinal);
%   O = make_i13ef_reviewer12_publication_outputs('/path/i13d_optimizer_benchmark.mat');
%
% Outputs include:
%   Fig_I13EF_optimizer_comparison_ieee.{pdf,png,fig}
%   Fig_I13EF_four_method_ieee.{pdf,png,fig}
%   Table_I13EF_optimizer_summary_by_M.csv
%   Table_I13EF_turbulent_method_comparison.csv
%   Table_I13EF_turbulent_gain_summary_by_M.csv
%   Table_I13EF_all_pointwise.csv
%   I13EF_reviewer12_summary.txt
%
% Scientific framing:
%   - Reviewer 1: generic literature GS (Zhou) vs channel-adaptive FSO GS.
%   - Reviewer 2: frozen SA vs constrained Pattern Search.
%   - No global-optimum claim is made.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'MakeFigures',true,@logical_scalar);
    addParameter(p,'CloseFigures',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    result=load_result(result);
    validate_result(result);

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=fullfile(result.outputDirectory,'publication_ready');
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    P=sortrows(result.physicalSummary,{'M','SNRdB','SigmaR2'});
    C=sortrows(result.caseSummary,{'M','SNRdB','SigmaR2','Replicate'});
    P.Outcome=string(P.Outcome);

    pointwise=build_pointwise_table(P);
    optimizerByM=build_optimizer_by_m(P,C);
    turbulent=pointwise(abs(pointwise.SigmaR2-0.3)<1e-12,:);
    turbulentGainByM=build_turbulent_gain_by_m(turbulent);

    paths=struct();
    paths.pointwise=fullfile(outDir,'Table_I13EF_all_pointwise.csv');
    paths.optimizerByM=fullfile(outDir,'Table_I13EF_optimizer_summary_by_M.csv');
    paths.turbulent=fullfile(outDir,'Table_I13EF_turbulent_method_comparison.csv');
    paths.turbulentGainByM=fullfile(outDir,'Table_I13EF_turbulent_gain_summary_by_M.csv');
    paths.summaryText=fullfile(outDir,'I13EF_reviewer12_summary.txt');

    writetable(pointwise,paths.pointwise);
    writetable(optimizerByM,paths.optimizerByM);
    writetable(turbulent,paths.turbulent);
    writetable(turbulentGainByM,paths.turbulentGainByM);

    summary=build_summary(result,P,C,pointwise,turbulent);
    write_summary_text(paths.summaryText,summary,optimizerByM,turbulentGainByM);

    figPaths=struct();
    if o.MakeFigures
        f1=make_optimizer_figure(P,result.config.tieToleranceBits);
        figPaths.optimizer=save_ieee_triplet(f1,outDir,'Fig_I13EF_optimizer_comparison_ieee');

        f2=make_four_method_figure(P);
        figPaths.fourMethod=save_ieee_triplet(f2,outDir,'Fig_I13EF_four_method_ieee');

        if o.CloseFigures
            close(f1);
            close(f2);
        end
    end

    out=struct();
    out.version='I13EF-v1';
    out.outputDirectory=outDir;
    out.pointwise=pointwise;
    out.optimizerSummaryByM=optimizerByM;
    out.turbulentComparison=turbulent;
    out.turbulentGainSummaryByM=turbulentGainByM;
    out.summary=summary;
    out.paths=paths;
    out.figures=figPaths;
    out.globalOptimumClaim=false;

    save(fullfile(outDir,'i13ef_publication_outputs.mat'),'out','-v7.3');

    print_summary(out);
end


function result=load_result(x)
    if ischar(x)||isstring(x)
        S=load(char(x));
        if isfield(S,'result')
            result=S.result;
        elseif isfield(S,'Rfinal')
            result=S.Rfinal;
        else
            error('make_i13ef_reviewer12_publication_outputs:MAT', ...
                'MAT file must contain result or Rfinal.');
        end
    else
        result=x;
    end
end


function validate_result(result)
    if ~isstruct(result) || ~isfield(result,'physicalSummary') || ...
            ~isfield(result,'caseSummary') || ~isfield(result,'config')
        error('make_i13ef_reviewer12_publication_outputs:Input', ...
            'Expected an I.13D result struct or MAT path.');
    end

    P=result.physicalSummary;
    needP={'M','SNRdB','SigmaR2','PAMValidated','ZhouValidated', ...
        'SAMean','PSMean','SAStd','PSStd','DeltaRep1','DeltaRep2', ...
        'SAMinusPSMean','Outcome'};
    missing=setdiff(needP,P.Properties.VariableNames);
    if ~isempty(missing)
        error('make_i13ef_reviewer12_publication_outputs:PhysicalColumns', ...
            'Missing physical-summary columns: %s',strjoin(missing,', '));
    end
    if any(~isfinite(P.ZhouValidated))
        error('make_i13ef_reviewer12_publication_outputs:ZhouMissing', ...
            'Zhou baseline values are required for the Reviewer-1/2 combined outputs.');
    end

    C=result.caseSummary;
    needC={'M','PSRestartStd','PSTotalFunctionEvaluations','PSMaxEvalOverrun'};
    missing=setdiff(needC,C.Properties.VariableNames);
    if ~isempty(missing)
        error('make_i13ef_reviewer12_publication_outputs:CaseColumns', ...
            'Missing case-summary columns: %s',strjoin(missing,', '));
    end
end


function T=build_pointwise_table(P)
    T=table();
    T.M=P.M;
    T.SNRdB=P.SNRdB;
    T.SigmaR2=P.SigmaR2;
    T.PAM=P.PAMValidated;
    T.ZhouGS=P.ZhouValidated;
    T.PatternSearchGS=P.PSMean;
    T.SAGS=P.SAMean;
    T.ZhouMinusPAM=T.ZhouGS-T.PAM;
    T.PatternSearchMinusZhou=T.PatternSearchGS-T.ZhouGS;
    T.SAMinusZhou=T.SAGS-T.ZhouGS;
    T.SAMinusPatternSearch=T.SAGS-T.PatternSearchGS;
    T.SAStdAcrossReplicates=P.SAStd;
    T.PSStdAcrossReplicates=P.PSStd;
    T.Outcome=string(P.Outcome);
end


function T=build_optimizer_by_m(P,C)
    MVec=unique(P.M).';
    n=numel(MVec);
    M=MVec(:);
    NPoints=zeros(n,1);
    Equivalent=zeros(n,1);
    SABetter=zeros(n,1);
    PatternSearchBetter=zeros(n,1);
    MeanSAMinusPS=zeros(n,1);
    MedianSAMinusPS=zeros(n,1);
    MinSAMinusPS=zeros(n,1);
    MaxSAMinusPS=zeros(n,1);
    MaxAbsSAMinusPS=zeros(n,1);
    MeanPSRestartStd=zeros(n,1);
    MaxPSRestartStd=zeros(n,1);
    MeanSAStdAcrossReplicates=zeros(n,1);
    MeanPSStdAcrossReplicates=zeros(n,1);

    for k=1:n
        p=P(P.M==M(k),:);
        c=C(C.M==M(k),:);
        outcome=string(p.Outcome);
        d=p.SAMinusPSMean;

        NPoints(k)=height(p);
        Equivalent(k)=sum(outcome=="Equivalent");
        SABetter(k)=sum(outcome=="SA better");
        PatternSearchBetter(k)=sum(outcome=="Pattern Search better");
        MeanSAMinusPS(k)=mean(d);
        MedianSAMinusPS(k)=median(d);
        MinSAMinusPS(k)=min(d);
        MaxSAMinusPS(k)=max(d);
        MaxAbsSAMinusPS(k)=max(abs(d));
        MeanPSRestartStd(k)=mean(c.PSRestartStd);
        MaxPSRestartStd(k)=max(c.PSRestartStd);
        MeanSAStdAcrossReplicates(k)=mean(p.SAStd);
        MeanPSStdAcrossReplicates(k)=mean(p.PSStd);
    end

    T=table(M,NPoints,Equivalent,SABetter,PatternSearchBetter, ...
        MeanSAMinusPS,MedianSAMinusPS,MinSAMinusPS,MaxSAMinusPS, ...
        MaxAbsSAMinusPS,MeanPSRestartStd,MaxPSRestartStd, ...
        MeanSAStdAcrossReplicates,MeanPSStdAcrossReplicates);
end


function T=build_turbulent_gain_by_m(A)
    MVec=unique(A.M).';
    n=numel(MVec);
    M=MVec(:);
    NPoints=zeros(n,1);
    MeanZhouMinusPAM=zeros(n,1);
    MeanPatternSearchMinusZhou=zeros(n,1);
    MeanSAMinusZhou=zeros(n,1);
    MeanSAMinusPatternSearch=zeros(n,1);

    for k=1:n
        q=A(A.M==M(k),:);
        NPoints(k)=height(q);
        MeanZhouMinusPAM(k)=mean(q.ZhouMinusPAM);
        MeanPatternSearchMinusZhou(k)=mean(q.PatternSearchMinusZhou);
        MeanSAMinusZhou(k)=mean(q.SAMinusZhou);
        MeanSAMinusPatternSearch(k)=mean(q.SAMinusPatternSearch);
    end

    T=table(M,NPoints,MeanZhouMinusPAM,MeanPatternSearchMinusZhou, ...
        MeanSAMinusZhou,MeanSAMinusPatternSearch);
end


function S=build_summary(result,P,C,pointwise,turbulent)
    d=P.SAMinusPSMean;
    outcome=string(P.Outcome);

    S=struct();
    S.nPhysicalPoints=height(P);
    S.nOptimizationCases=height(C);
    S.nEquivalent=sum(outcome=="Equivalent");
    S.nSABetter=sum(outcome=="SA better");
    S.nPatternSearchBetter=sum(outcome=="Pattern Search better");
    S.meanSAMinusPS=mean(d);
    S.medianSAMinusPS=median(d);
    S.minSAMinusPS=min(d);
    S.maxSAMinusPS=max(d);
    S.maxAbsSAMinusPS=max(abs(d));

    p0=P(abs(P.SigmaR2)<1e-12,:);
    pT=P(abs(P.SigmaR2-0.3)<1e-12,:);
    S.meanDeltaSigma0=mean(p0.SAMinusPSMean);
    S.meanDeltaSigma03=mean(pT.SAMinusPSMean);

    S.meanTurbulentZhouGain=mean(turbulent.ZhouMinusPAM);
    S.meanTurbulentPSAdaptationGain=mean(turbulent.PatternSearchMinusZhou);
    S.meanTurbulentSAAdaptationGain=mean(turbulent.SAMinusZhou);
    S.meanTurbulentSAMinusPS=mean(turbulent.SAMinusPatternSearch);

    if isfield(result,'metrics') && isfield(result.metrics,'totalPatternSearchFunctionEvaluations')
        S.actualPSSearchEvaluations=result.metrics.totalPatternSearchFunctionEvaluations;
    else
        S.actualPSSearchEvaluations=sum(C.PSTotalFunctionEvaluations);
    end

    if isfield(result.config,'totalMaxFunctionEvaluations')
        S.maximumPSSearchEvaluationBudget=result.config.totalMaxFunctionEvaluations;
    else
        S.maximumPSSearchEvaluationBudget=NaN;
    end

    if isfinite(S.maximumPSSearchEvaluationBudget)
        S.psBudgetUtilization=S.actualPSSearchEvaluations/S.maximumPSSearchEvaluationBudget;
    else
        S.psBudgetUtilization=NaN;
    end

    S.maxPSEvaluationOverrun=max(C.PSMaxEvalOverrun);
    S.tieToleranceBits=result.config.tieToleranceBits;

    [~,iMin]=min(d);
    [~,iMax]=max(d);
    S.patternSearchBestPoint=point_descriptor(P(iMin,:));
    S.patternSearchBestDelta=d(iMin);
    S.saBestPoint=point_descriptor(P(iMax,:));
    S.saBestDelta=d(iMax);
end


function s=point_descriptor(R)
    s=sprintf('M=%d, SNR=%g dB, sigma_R^2=%.1f',R.M,R.SNRdB,R.SigmaR2);
end


function write_summary_text(path,S,optimizerByM,turbulentGainByM)
    fid=fopen(path,'w');
    if fid<0,error('make_i13ef_reviewer12_publication_outputs:File','Cannot open %s.',path);end
    cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid,'I.13E-F Reviewer 1/2 publication summary\n');
    fprintf(fid,'=========================================\n\n');
    fprintf(fid,'No global-optimum claim is made. All AMI values are validated results.\n\n');

    fprintf(fid,'Reviewer 2 -- optimizer robustness\n');
    fprintf(fid,'Physical points: %d | optimization cases: %d\n',S.nPhysicalPoints,S.nOptimizationCases);
    fprintf(fid,'Predeclared equivalence tolerance: %.1e bit/symbol\n',S.tieToleranceBits);
    fprintf(fid,'Equivalent: %d/%d | SA better: %d/%d | Pattern Search better: %d/%d\n', ...
        S.nEquivalent,S.nPhysicalPoints,S.nSABetter,S.nPhysicalPoints, ...
        S.nPatternSearchBetter,S.nPhysicalPoints);
    fprintf(fid,'SA-PS mean / median: %+.6f / %+.6f bit/symbol\n', ...
        S.meanSAMinusPS,S.medianSAMinusPS);
    fprintf(fid,'SA-PS min / max: %+.6f / %+.6f bit/symbol\n', ...
        S.minSAMinusPS,S.maxSAMinusPS);
    fprintf(fid,'Mean SA-PS at sigma_R^2=0: %.6f bit/symbol\n',S.meanDeltaSigma0);
    fprintf(fid,'Mean SA-PS at sigma_R^2=0.3: %.6f bit/symbol\n',S.meanDeltaSigma03);
    fprintf(fid,'Pattern Search best relative point: %s | SA-PS=%+.6f\n', ...
        S.patternSearchBestPoint,S.patternSearchBestDelta);
    fprintf(fid,'SA largest relative advantage: %s | SA-PS=%+.6f\n', ...
        S.saBestPoint,S.saBestDelta);
    fprintf(fid,'Pattern Search search evaluations: %.0f',S.actualPSSearchEvaluations);
    if isfinite(S.maximumPSSearchEvaluationBudget)
        fprintf(fid,' / %.0f maximum (%.1f%% used)\n', ...
            S.maximumPSSearchEvaluationBudget,100*S.psBudgetUtilization);
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'Maximum PS evaluation-budget overrun: %+g\n\n',S.maxPSEvaluationOverrun);

    fprintf(fid,'Reviewer 1 -- literature GS vs channel-adaptive FSO GS, sigma_R^2=0.3\n');
    fprintf(fid,'Mean Zhou-PAM generic GS gain: %.6f bit/symbol\n',S.meanTurbulentZhouGain);
    fprintf(fid,'Mean PatternSearch-Zhou adaptation gain: %.6f bit/symbol\n', ...
        S.meanTurbulentPSAdaptationGain);
    fprintf(fid,'Mean SA-Zhou adaptation gain: %.6f bit/symbol\n', ...
        S.meanTurbulentSAAdaptationGain);
    fprintf(fid,'Mean SA-PatternSearch gap: %.6f bit/symbol\n\n', ...
        S.meanTurbulentSAMinusPS);

    fprintf(fid,'Optimizer summary by M\n');
    fprintf(fid,'M  N  Eq  SA>PS  PS>SA  mean(SA-PS)  median(SA-PS)\n');
    for k=1:height(optimizerByM)
        r=optimizerByM(k,:);
        fprintf(fid,'%d  %d  %d  %d  %d  %+.6f  %+.6f\n', ...
            r.M,r.NPoints,r.Equivalent,r.SABetter,r.PatternSearchBetter, ...
            r.MeanSAMinusPS,r.MedianSAMinusPS);
    end

    fprintf(fid,'\nTurbulent gain summary by M (sigma_R^2=0.3)\n');
    fprintf(fid,'M  Zhou-PAM  PS-Zhou  SA-Zhou  SA-PS\n');
    for k=1:height(turbulentGainByM)
        r=turbulentGainByM(k,:);
        fprintf(fid,'%d  %+.6f  %+.6f  %+.6f  %+.6f\n', ...
            r.M,r.MeanZhouMinusPAM,r.MeanPatternSearchMinusZhou, ...
            r.MeanSAMinusZhou,r.MeanSAMinusPatternSearch);
    end
end


function fig=make_optimizer_figure(P,tol)
    P=sortrows(P,{'M','SNRdB','SigmaR2'});
    MVec=unique(P.M).';
    colors=lines(numel(MVec));

    fig=figure('Color','w','Name','I13EF optimizer comparison', ...
        'Units','inches','Position',[1 1 7.16 3.15]);
    tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    % Panel (a): direct SA-vs-PS agreement.
    ax=nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    allVals=[P.SAMean;P.PSMean];
    lo=min(allVals); hi=max(allVals);
    pad=max(0.03,0.035*(hi-lo+eps));
    hDiag=plot(ax,[lo-pad hi+pad],[lo-pad hi+pad],'k--','LineWidth',1.0);

    hM=gobjects(numel(MVec),1);
    for iM=1:numel(MVec)
        rows=find(P.M==MVec(iM));
        for q=1:numel(rows)
            k=rows(q);
            if abs(P.SigmaR2(k))<1e-12,marker='o';else,marker='s';end
            if P.SNRdB(k)>=30,face=colors(iM,:);else,face='w';end
            h=scatter(ax,P.SAMean(k),P.PSMean(k),38, ...
                'Marker',marker,'MarkerEdgeColor',colors(iM,:), ...
                'MarkerFaceColor',face,'LineWidth',1.0);
            if q==1,hM(iM)=h;end
        end
    end
    xlabel(ax,'SA validated AMI [bits/symbol]');
    ylabel(ax,'Pattern Search validated AMI [bits/symbol]');
    title(ax,'(a) Optimizer agreement','FontWeight','normal');
    xlim(ax,[lo-pad hi+pad]); ylim(ax,[lo-pad hi+pad]);
    axis(ax,'square');
    legLabels=arrayfun(@(m)sprintf('M=%d',m),MVec,'UniformOutput',false);
    legend(ax,[hDiag;hM],[{'y=x'},legLabels],'Location','southeast','FontSize',7);
    text(ax,0.03,0.97,'o: sigma_R^2=0;  s: sigma_R^2=0.3;  filled: 30 dB', ...
        'Units','normalized','VerticalAlignment','top','FontSize',6.6, ...
        'Interpreter','none');
    apply_ieee_axes(ax);

    % Panel (b): delta grouped by M and channel condition.
    ax=nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    xM=1:numel(MVec);
    bandX=[0.55 numel(MVec)+0.45 numel(MVec)+0.45 0.55];
    patch(ax,bandX,[-tol -tol tol tol],[0.90 0.90 0.90], ...
        'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    yline(ax,0,'k-','LineWidth',0.8,'HandleVisibility','off');
    yline(ax,tol,'k:','LineWidth',0.7,'HandleVisibility','off');
    yline(ax,-tol,'k:','LineWidth',0.7,'HandleVisibility','off');

    cond=[20 0;20 0.3;30 0;30 0.3];
    offsets=[-0.21 -0.07 0.07 0.21];
    markers={'o','s','o','s'};
    fills=[false false true true];
    ccol=lines(4);
    hCond=gobjects(4,1);

    for j=1:4
        xx=nan(numel(MVec),1); yy=xx; ee=xx;
        for iM=1:numel(MVec)
            idx=P.M==MVec(iM) & P.SNRdB==cond(j,1) & abs(P.SigmaR2-cond(j,2))<1e-12;
            r=P(idx,:);
            xx(iM)=iM+offsets(j);
            yy(iM)=r.SAMinusPSMean;
            ee(iM)=std([r.DeltaRep1 r.DeltaRep2],0,2);
        end
        if fills(j),face=ccol(j,:);else,face='w';end
        h=errorbar(ax,xx,yy,ee,'LineStyle','none','Color',ccol(j,:), ...
            'Marker',markers{j},'MarkerSize',5.0,'MarkerFaceColor',face, ...
            'MarkerEdgeColor',ccol(j,:),'LineWidth',0.9,'CapSize',3);
        hCond(j)=h;
    end

    xlim(ax,[0.5 numel(MVec)+0.5]);
    xticks(ax,xM); xticklabels(ax,arrayfun(@num2str,MVec,'UniformOutput',false));
    xlabel(ax,'Constellation order M');
    ylabel(ax,'I_{SA}-I_{PS} [bits/symbol]');
    title(ax,'(b) Difference by operating condition','FontWeight','normal');
    labels={'20 dB, sigma_R^2=0','20 dB, sigma_R^2=0.3', ...
        '30 dB, sigma_R^2=0','30 dB, sigma_R^2=0.3'};
    legend(ax,hCond,labels,'Location','northwest','FontSize',6.7);
    apply_ieee_axes(ax);
end


function fig=make_four_method_figure(P)
    P=sortrows(P,{'M','SNRdB','SigmaR2'});
    MVec=unique(P.M).';
    cond=[20 0;20 0.3;30 0;30 0.3];
    x=1:4;
    offsets=[-0.18 -0.06 0.06 0.18];
    methodColors=lines(4);
    markers={'o','d','s','^'};
    labels={'Uniform PAM','Zhou GS','Pattern Search','SA'};

    fig=figure('Color','w','Name','I13EF four-method comparison', ...
        'Units','inches','Position',[1 1 7.16 5.25]);
    tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    hLegend=gobjects(4,1);
    for iM=1:numel(MVec)
        R=P(P.M==MVec(iM),:);
        ax=nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

        pam=nan(4,1); zhou=pam; ps=pam; sa=pam; psErr=pam; saErr=pam;
        for j=1:4
            idx=R.SNRdB==cond(j,1) & abs(R.SigmaR2-cond(j,2))<1e-12;
            q=R(idx,:);
            pam(j)=q.PAMValidated;
            zhou(j)=q.ZhouValidated;
            ps(j)=q.PSMean;
            sa(j)=q.SAMean;
            psErr(j)=q.PSStd;
            saErr(j)=q.SAStd;
        end

        vals={pam,zhou,ps,sa};
        errs={zeros(4,1),zeros(4,1),psErr,saErr};
        for m=1:4
            xx=x+offsets(m);
            if m<=2
                h=plot(ax,xx,vals{m},'LineStyle','none','Marker',markers{m}, ...
                    'MarkerSize',5.0,'MarkerEdgeColor',methodColors(m,:), ...
                    'MarkerFaceColor','w','LineWidth',1.0);
            else
                h=errorbar(ax,xx,vals{m},errs{m},'LineStyle','none', ...
                    'Color',methodColors(m,:),'Marker',markers{m}, ...
                    'MarkerSize',5.0,'MarkerEdgeColor',methodColors(m,:), ...
                    'MarkerFaceColor',methodColors(m,:),'LineWidth',0.9,'CapSize',3);
            end
            if iM==1,hLegend(m)=h;end
        end

        xlim(ax,[0.55 4.45]);
        xticks(ax,x);
        xticklabels(ax,{'20 / 0','20 / 0.3','30 / 0','30 / 0.3'});
        xlabel(ax,'SNR [dB] / \\sigma_R^2','Interpreter','tex');
        ylabel(ax,'Validated AMI [bits/symbol]');
        title(ax,sprintf('M = %d',MVec(iM)),'FontWeight','normal');
        apply_ieee_axes(ax);
    end

    lg=legend(hLegend,labels,'Orientation','horizontal','FontSize',7);
    lg.Layout.Tile='north';
end


function apply_ieee_axes(ax)
    ax.Color='w';
    ax.XColor='k';
    ax.YColor='k';
    ax.GridColor=[0.82 0.82 0.82];
    ax.GridAlpha=0.55;
    ax.LineWidth=0.75;
    ax.FontName='Times New Roman';
    ax.FontSize=8;
    try,ax.Toolbar.Visible='off';catch,end
end


function p=save_ieee_triplet(fig,outDir,stem)
    p=struct();
    p.png=fullfile(outDir,[stem '.png']);
    p.pdf=fullfile(outDir,[stem '.pdf']);
    p.fig=fullfile(outDir,[stem '.fig']);

    if ~isgraphics(fig,'figure')
        error('make_i13ef_reviewer12_publication_outputs:Figure','Invalid figure handle.');
    end

    drawnow;
    try
        savefig(fig,p.fig);
    catch
        hgsave(fig,p.fig);
    end

    exportgraphics(fig,p.png,'Resolution',400,'BackgroundColor','white');
    try
        exportgraphics(fig,p.pdf,'ContentType','vector','BackgroundColor','white');
    catch
        exportgraphics(fig,p.pdf,'BackgroundColor','white');
    end
end


function print_summary(out)
    S=out.summary;
    fprintf('\n============================================================\n');
    fprintf('I.13E-F REVIEWER 1/2 PUBLICATION OUTPUTS COMPLETE\n');
    fprintf('Equivalent: %d/%d | SA better: %d/%d | PS better: %d/%d\n', ...
        S.nEquivalent,S.nPhysicalPoints,S.nSABetter,S.nPhysicalPoints, ...
        S.nPatternSearchBetter,S.nPhysicalPoints);
    fprintf('Mean SA-PS: %+.6f bit/symbol | sigma_R^2=0: %.6f | sigma_R^2=0.3: %.6f\n', ...
        S.meanSAMinusPS,S.meanDeltaSigma0,S.meanDeltaSigma03);
    fprintf('Turbulent mean gains: Zhou-PAM=%.6f | PS-Zhou=%.6f | SA-Zhou=%.6f\n', ...
        S.meanTurbulentZhouGain,S.meanTurbulentPSAdaptationGain,S.meanTurbulentSAAdaptationGain);
    fprintf('Output: %s\n',out.outputDirectory);
    fprintf('No optimization was rerun. No global-optimum claim is made.\n');
    fprintf('============================================================\n');
end


function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
