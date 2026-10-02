function out = make_ieee_revision_figures(varargin)
%MAKE_IEEE_REVISION_FIGURES I.12B: regenerate manuscript Figs. 3-5.
% Reads the frozen I.11.1 production bundle only; no SA/AMI recomputation.
%
% Performance curves use the mean validated AMI across two replicates.
% Fig. 4 uses replicate 1 consistently for geometry (no best-replicate pick).

p=inputParser; p.FunctionName=mfilename;
addParameter(p,'ProductionBundlePath','',@istext);
addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@istext);
addParameter(p,'OutputDirectory','',@istext);
addParameter(p,'ShowReplicateErrorBars',false,@islogicalscalar);
parse(p,varargin{:}); o=p.Results;

A=analyze_i11p1_production_results( ...
    'ProductionBundlePath',o.ProductionBundlePath, ...
    'ResultsRoot',o.ResultsRoot,'WriteCSV',true);

if strlength(string(o.OutputDirectory))>0
    outDir=char(string(o.OutputDirectory));
else
    outDir=fullfile(A.outputDirectory,'ieee_figures');
end
if exist(outDir,'dir')~=7,mkdir(outDir);end
P=A.physicalTable;
assert(height(P)==96,'Expected 96 physical points.');

f3=fig3(P,o.ShowReplicateErrorBars);
f4=fig4(A.fig4_M8_SNR20);
f5=fig5(P,o.ShowReplicateErrorBars);

out=struct();
out.version='I12B-v1';
out.outputDirectory=outDir;
out.fig3=saveall(f3,outDir,'Fig3_M8_AMI_vs_SNR');
out.fig4=saveall(f4,outDir,'Fig4_M8_SNR20_constellations');
out.fig5=saveall(f5,outDir,'Fig5_AMI_vs_SNR_allM');
out.tableI=fullfile(outDir,'TableI_M8_SNR20_rows.tex');
write_table(A.tableI_M8_SNR20,out.tableI);
out.digest=fullfile(outDir,'I12B_results_digest.txt');
write_digest(A,out.digest);
save(fullfile(outDir,'i12b_figure_manifest.mat'),'out','-v7.3');

fprintf('\nI.12B complete. Outputs: %s\n',outDir);
fprintf('Fig3 PDF: %s\n',out.fig3.pdf);
fprintf('Fig4 PDF: %s\n',out.fig4.pdf);
fprintf('Fig5 PDF: %s\n',out.fig5.pdf);
end

function f=fig3(P,showErr)
M=8; sig=unique(P.SigmaX2(P.M==M)).'; C=lines(numel(sig));
f=figure('Color','w','Name','Fig3','Units','inches','Position',[1 1 7 4.6]);
ax=axes(f); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for i=1:numel(sig)
    R=sortrows(P(P.M==M & abs(P.SigmaX2-sig(i))<1e-12,:),'SNRdB');
    plot(ax,R.SNRdB,R.PAMValidated,'--','Color',C(i,:),'LineWidth',1.2, ...
        'DisplayName',sprintf('PAM, \\sigma_X^2=%.1f',sig(i)));
    if showErr
        errorbar(ax,R.SNRdB,R.GSMeanValidated,R.GSStdAcrossReplicates,'-o', ...
            'Color',C(i,:),'LineWidth',1.6,'MarkerSize',4,'CapSize',3, ...
            'DisplayName',sprintf('GS, \\sigma_X^2=%.1f',sig(i)));
    else
        plot(ax,R.SNRdB,R.GSMeanValidated,'-o','Color',C(i,:), ...
            'LineWidth',1.6,'MarkerSize',4, ...
            'DisplayName',sprintf('GS, \\sigma_X^2=%.1f',sig(i)));
    end
end
yline(ax,log2(M),':','HandleVisibility','off');
xlabel(ax,'SNR [dB]'); ylabel(ax,'AMI [bits/symbol]');
xlim(ax,[5 30]); ylim(ax,[0 3.1]);
legend(ax,'Location','best','NumColumns',2,'FontSize',8);
set(ax,'FontSize',9);
end

function f=fig4(T)
T=sortrows(T,'SigmaX2'); M=8; k=1:M;
f=figure('Color','w','Name','Fig4','Units','inches','Position',[1 1 7 5.2]);
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for i=1:height(T)
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    pam=zeros(1,M); gs=zeros(1,M);
    for j=1:M
        pam(j)=T.(sprintf('PAM_x%d',j))(i);
        gs(j)=T.(sprintf('GS_representative_x%d',j))(i);
    end
    plot(ax,k,pam,'--o','LineWidth',1.2,'MarkerSize',4,'DisplayName','Uniform PAM');
    plot(ax,k,gs,'-s','LineWidth',1.6,'MarkerSize',4,'DisplayName','GS (rep. 1)');
    title(ax,sprintf('\\sigma_X^2 = %.1f',T.SigmaX2(i)),'FontWeight','normal');
    xlabel(ax,'Symbol index'); ylabel(ax,'Intensity level'); xlim(ax,[1 M]); xticks(ax,1:M);
    if i==1,legend(ax,'Location','best','FontSize',8);end
    set(ax,'FontSize',9);
end
end

function f=fig5(P,showErr)
Mvec=unique(P.M).'; sig=unique(P.SigmaX2).'; C=lines(numel(sig));
f=figure('Color','w','Name','Fig5','Units','inches','Position',[1 1 7.1 6.2]);
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for im=1:numel(Mvec)
    M=Mvec(im); ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for i=1:numel(sig)
        R=sortrows(P(P.M==M & abs(P.SigmaX2-sig(i))<1e-12,:),'SNRdB');
        hv=onoff(im==1);
        plot(ax,R.SNRdB,R.PAMValidated,'--','Color',C(i,:),'LineWidth',1, ...
            'HandleVisibility',hv,'DisplayName',sprintf('PAM, \\sigma_X^2=%.1f',sig(i)));
        if showErr
            errorbar(ax,R.SNRdB,R.GSMeanValidated,R.GSStdAcrossReplicates,'-o', ...
                'Color',C(i,:),'LineWidth',1.3,'MarkerSize',3,'CapSize',2, ...
                'HandleVisibility',hv,'DisplayName',sprintf('GS, \\sigma_X^2=%.1f',sig(i)));
        else
            plot(ax,R.SNRdB,R.GSMeanValidated,'-o','Color',C(i,:),'LineWidth',1.3, ...
                'MarkerSize',3,'HandleVisibility',hv, ...
                'DisplayName',sprintf('GS, \\sigma_X^2=%.1f',sig(i)));
        end
    end
    yline(ax,log2(M),':','HandleVisibility','off');
    title(ax,sprintf('M = %d',M),'FontWeight','normal');
    xlabel(ax,'SNR [dB]'); ylabel(ax,'AMI [bits/symbol]');
    xlim(ax,[5 30]); ylim(ax,[0 log2(M)*1.04]); set(ax,'FontSize',8.5);
    if im==1,legend(ax,'Location','best','NumColumns',2,'FontSize',7);end
end
end

function p=saveall(f,d,stem)
p.png=fullfile(d,[stem '.png']); p.pdf=fullfile(d,[stem '.pdf']); p.fig=fullfile(d,[stem '.fig']);
exportgraphics(f,p.png,'Resolution',300);
try,exportgraphics(f,p.pdf,'ContentType','vector');catch,exportgraphics(f,p.pdf);end
savefig(f,p.fig);
end

function write_table(T,path)
T=sortrows(T,'SigmaX2'); fid=fopen(path,'w'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% GS = mean validated AMI across two replicates.\n');
for i=1:height(T)
    fprintf(fid,'%.2f & %.4f & %.4f & %.4f \\\\\n', ...
        T.SigmaX2(i),T.PAMValidated(i),T.GSMeanValidated(i),T.GainMeanBits(i));
end
end

function write_digest(A,path)
M=A.metrics; T=A.tableI_M8_SNR20;
fid=fopen(path,'w'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'I.12B production digest\nSource: %s\n',A.sourceBundlePath);
fprintf(fid,'Selection changes: %d/%d (%.2f%%%%)\n',M.nSelectionChanges,M.nOptimizationCases,100*M.selectionChangeFraction);
fprintf(fid,'Max fast-validator gap: %.6e bits/symbol\n',M.maxFastValidatorGap);
fprintf(fid,'Max replicate AMI range: %.6f bits/symbol\n',M.maxReplicateRangeBits);
fprintf(fid,'Material negative gains (< -%.1e): %d individual, %d physical means\n\n', ...
    M.numericalZeroTolerance,M.nMateriallyNegativeIndividualGains,M.nMateriallyNegativeMeanPhysicalGains);
fprintf(fid,'M=8, SNR=20 dB: sigma_X^2, PAM, GS mean, gain, GS std\n');
for i=1:height(T)
    fprintf(fid,'%.2f, %.6f, %.6f, %.6f, %.6e\n',T.SigmaX2(i),T.PAMValidated(i), ...
        T.GSMeanValidated(i),T.GainMeanBits(i),T.GSStdAcrossReplicates(i));
end
end

function s=onoff(tf),if tf,s='on';else,s='off';end,end
function tf=istext(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=islogicalscalar(v),tf=islogical(v)&&isscalar(v);end
