function K52 = make_k5p2_fig7(varargin)
%MAKE_K5P2_FIG7  Commit K5.2: reproducible 0.3-turbulence Figure 7.
%
% DEFAULT (NO optimization, NO AMI evaluation):
%   S = make_k5p2_fig7('DryRun',true);
%   K52 = make_k5p2_fig7();
%
% OPTIONAL, only for independent Pattern Search re-check (12 tasks):
%   K52 = make_k5p2_fig7('RerunPatternSearch',true);
%
% Frozen I.11.1, I.13B, I.13D and K5 scientific data are READ-ONLY.
% Four-method primary plot always uses the original validated I.13D
% values, even when the optional PS re-check is requested.
%
% All output goes under results/ieee_revision2/commit_k5p2_fig7/.
% No global-optimum claim. This is a single representative test point.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@is_text);
    addParameter(p,'ZhouResultPath','',@is_text);
    addParameter(p,'I13DResultPath','',@is_text);
    addParameter(p,'K5ResultPath','',@is_text);
    addParameter(p,'OutputDirectory','',@is_text);
    addParameter(p,'ResumeDirectory','',@is_text);
    addParameter(p,'RerunPatternSearch',false,@is_bool);
    addParameter(p,'ParallelWorkers',8,@is_positive_int);
    addParameter(p,'DryRun',false,@is_bool);
    parse(p,varargin{:}); o=p.Results;
    if strlength(string(o.ResumeDirectory))>0 && ...
            strlength(string(o.OutputDirectory))>0
        error('K5p2:Directory','Choose OutputDirectory OR ResumeDirectory.');
    end

    [prodPath,zhouPath,i13dPath,k5Path]=resolve_sources(o);
    S=load(prodPath,'bundle'); B=S.bundle;
    if ~isfield(B,'productionFreeze') || ...
            ~strcmp(B.productionFreeze.freezeVersion,'I11.1-v1')
        error('K5p2:Freeze','Expected frozen I.11.1 source.');
    end
    Z=load(zhouPath,'result'); Z=Z.result;
    D=load(i13dPath,'result'); R=D.result;
    H=load(k5Path,'R5'); K5=H.R5;
    assert(strcmp(Z.version,'I13B-v1') && strcmp(R.version,'I13D-v1'), ...
        'K5p2:SourceVersions','Expected frozen I13B-v1/I13D-v1.');
    assert(same_path(prodPath,Z.sourceProduction) && ...
        same_path(prodPath,R.sourceProduction) && ...
        same_path(prodPath,K5.sourceProduction) && ...
        same_path(zhouPath,K5.sourceZhou), ...
        'K5p2:Provenance','All sources must refer to SAME I.11.1 production.');

    F=i13d_benchmark_config(); C=i11p1_production_config();
    verify_frozen(F,C,R.config);
    F.MVec=8;F.SNRVec=20;F.TurbulenceVec=0.3;
    F.RestartsByM=6;F.MaxIterByM=4000;F.MaxEvalsPerStart=4001;
    F.nPhysicalPoints=1;F.nOptimizationCases=2;
    F.nRestartTasks=12;F.totalMaxFunctionEvaluations=48012;
    [casePlan,taskPlan]=i13d_build_task_plan(B.summary.caseTable,F);
    if height(casePlan)~=2 || height(taskPlan)~=12 || ...
        ~isequal(sort(casePlan.Replicate).',[1 2]) || ...
        any(casePlan.NStarts~=6 | casePlan.MaxIterSA~=4000 | ...
            casePlan.MaxEvalsPerStart~=4001 | ...
            abs(casePlan.SelectedDt-0.01)>1e-12)
        error('K5p2:CasePlan','Frozen 0.3 task plan does not match K5 protocol.');
    end
    P=R.physicalSummary;
    q=P.M==8 & P.SNRdB==20 & abs(P.SigmaR2-0.3)<1e-12;
    if nnz(q)~=1,error('K5p2:PhysicalPoint','Exactly one M8/20/0.3 is required.');end
    point=P(q,:);

    origCases=R.caseSummary;
    idx=origCases.M==8 & origCases.SNRdB==20 & ...
        abs(origCases.SigmaR2-0.3)<1e-12;
    cases=sortrows(origCases(idx,:),'Replicate');
    if height(cases)~=2 || ~isequal(cases.Replicate(:).',[1 2])
        error('K5p2:Replicate','Expected two saved I.13D replicates.');
    end
    expected={'Seed','NStarts','SelectedDt','SAValidated', ...
        'PSValidatedCase','PAMValidated','PSBestRestart','MaxEvalsPerStart'};
    if ~isempty(setdiff(expected,cases.Properties.VariableNames))
        error('K5p2:CaseSchema','Missing I13D source case columns.');
    end
    if any(cases.Seed~=casePlan.Seed) || ...
            any(cases.NStarts~=6) || any(cases.MaxEvalsPerStart~=4001) || ...
            any(abs(cases.SelectedDt-0.01)>1e-12) || ...
            any(abs(cases.SAValidated-casePlan.SAValidated)>1e-9) || ...
            any(abs(cases.PAMValidated-casePlan.PAMValidated)>1e-9)
        error('K5p2:FrozenMismatch','Replicate seeds, AMI or budget mismatch.');
    end

    restarts=R.restartSummary;
    mask=restarts.M==8 & restarts.SNRdB==20 & ...
        abs(restarts.SigmaR2-0.3)<1e-12;
    restarts=sortrows(restarts(mask,:),{'Replicate','Restart'});
    if height(restarts)~=12 || numel(unique(restarts.IdentityKey))~=12 || ...
            any(restarts.EvalBudget~=4001) || ...
            any(restarts.FunctionEvaluations>4001) || ...
            any(~isfinite(restarts.PSValidated)) || ...
            any(abs(restarts.FastValidatorGap)>1e-3) || ...
            ~isequal(sort(string(restarts.IdentityKey)), ...
                     sort(string(taskPlan.IdentityKey)))
        error('K5p2:RestartAudit','I13D restart set or validation/budget failed.');
    end
    psByRep=zeros(2,1);
    for rep=1:2
        row=restarts(restarts.Replicate==rep,:);
        if height(row)~=6 || ~isequal(sort(row.Restart).',1:6)
            error('K5p2:RestartCount','Incomplete original restart group.');
        end
        [psByRep(rep),ii]=max(row.PSValidated);
        if abs(psByRep(rep)-cases.PSValidatedCase(rep))>1e-9 || ...
                row.Restart(ii)~=cases.PSBestRestart(rep)
            error('K5p2:Winner','Original I13D validated winner inconsistent.');
        end
    end
    pam=mean(cases.PAMValidated);
    sa=mean(cases.SAValidated);
    ps=mean(psByRep);
    saStd=std(cases.SAValidated,0);
    psStd=std(psByRep,0);

    z=Z.comparison;
    zi=z.M==8 & z.SNRdB==20 & abs(z.SigmaX2-0.3)<1e-12;
    if nnz(zi)~=1,error('K5p2:Zhou','Expected one I.13B comparison row.');end
    zhou=z.ZhouValidated(zi);

    numeric=[point.PAMValidated,point.SAMean,point.PSMean, ...
             point.SAStd,point.PSStd,point.SARep1,point.SARep2, ...
             point.PSRep1,point.PSRep2,point.ZhouValidated];
    expectedNumeric=[pam,sa,ps,saStd,psStd,cases.SAValidated.',psByRep.',zhou];
    if any(~isfinite(numeric)) || ...
            any(abs(numeric-expectedNumeric)>1e-9) || ...
            abs(z.PAMValidated(zi)-pam)>1e-9 || ...
            abs(z.GSMeanValidated(zi)-sa)>1e-9
        error('K5p2:Aggregation','Frozen AMI does not reproduce from replicates.');
    end
    if ps>=sa
        fprintf('K5p2 NOTE: SA is not ahead at this point; no cherry-picked assertion.\n');
    end

    methods=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    value=[pam;zhou;ps;sa];
    sampleStd=[0;0;psStd;saStd];
    four=table(methods,value,sampleStd, ...
        'VariableNames',{'Method','ValidatedAMI','StdAcrossReplicates'});

    if ~isfield(K5,'fourMethodTable') || ...
            ~isequal(string(K5.fourMethodTable.Method(:)),methods) || ...
            numel(K5.replicatePatternSearchAMI)~=2 || ...
            any(K5.casePlan.SigmaR2~=0.1) || ...
            numel(unique(K5.restartSummary.Restart))~=6 || ...
            height(K5.restartSummary)~=12
        error('K5p2:K5','Incomplete/incorrect previous K5 output.');
    end
    prev=K5.fourMethodTable;
    oldVals=prev.ValidatedAMI;
    z01=z.M==8 & z.SNRdB==20 & abs(z.SigmaX2-0.1)<1e-12;
    if nnz(z01)~=1 || ...
            abs(oldVals(1)-mean(K5.casePlan.PAMValidated))>1e-9 || ...
            abs(oldVals(2)-z.ZhouValidated(z01))>1e-9 || ...
            abs(oldVals(4)-mean(K5.casePlan.SAValidated))>1e-9 || ...
            abs(oldVals(3)-mean(K5.replicatePatternSearchAMI))>1e-9 || ...
            abs(prev.StdAcrossReplicates(3)-std(K5.replicatePatternSearchAMI,0))>1e-9
        error('K5p2:OldK5','Previous K5 replicate aggregation is inconsistent.');
    end
    for rep=1:2
        k5r=K5.restartSummary(K5.restartSummary.Replicate==rep,:);
        if height(k5r)~=6 || ...
                abs(max(k5r.PSValidated)-K5.replicatePatternSearchAMI(rep))>1e-9
            error('K5p2:OldK5Restarts','K5 validated restart winner mismatch.');
        end
    end
    comparison=table(methods,oldVals,value,value-oldVals, ...
        'VariableNames',{'Method','K5_sigmaR2_0p1', ...
        'K5p2_sigmaR2_0p3','Difference_0p3_minus_0p1'});
    gOld=oldVals(4)-oldVals(3);
    gNew=value(4)-value(3);
    gap=table([0.1;0.3],[gOld;gNew], ...
        'VariableNames',{'SigmaR2','SAMinusPatternSearch'});
    ratio=gNew/gOld;
    if ~isfinite(ratio) || gOld<=0
        error('K5p2:Gap','Cannot compute a positive gap ratio.');
    end

    K52=struct('version','K5.2-v1','dryRun',logical(o.DryRun), ...
        'M',8,'SNRdB',20,'SigmaR2',0.3,'sourceProduction',prodPath, ...
        'sourceZhou',zhouPath,'sourceI13D',i13dPath,'sourceK5',k5Path, ...
        'config',F,'casePlan',casePlan,'taskPlan',taskPlan, ...
        'frozenRestartRows',restarts,'frozenCaseRows',cases, ...
        'frozenPhysicalRow',point,'fourMethodTable',four, ...
        'k5Comparison',comparison,'optimizerGapByTurbulence',gap, ...
        'gapRatio',ratio,'globalOptimumClaim',false, ...
        'optimizationRerun',false,'figureSource','frozen I.13D');

    fprintf('\nCOMMIT K5.2 | M=8 SNR=20 sigma_R^2=0.3\n');
    fprintf('Source: frozen I.13D, 2 replicates, 12 validated restarts.\n');
    fprintf('SA-PS: %.8f at 0.1 -> %.8f at 0.3 | ratio %.2fx\n',gOld,gNew,ratio);
    disp(four);

    if o.DryRun
        fprintf('DRY RUN: no output written; no optimization, no AMI evaluation.\n');
        return;
    end
    outDir=prepare_directory(o);
    K52.outputDirectory=outDir;
    writetable(four,fullfile(outDir,'Fig7_K5p2_four_methods.csv'));
    writetable(restarts,fullfile(outDir,'Fig7_K5p2_frozen_restarts.csv'));
    writetable(cases,fullfile(outDir,'Fig7_K5p2_replicates.csv'));
    writetable(comparison,fullfile(outDir,'K5_vs_K5p2_comparison.csv'));
    writetable(gap,fullfile(outDir,'K5_vs_K5p2_SA_PS_gap.csv'));
    K52.figurePaths=plot_k5p2_fig7_four_methods(four,0.3,outDir,'CloseFigure',true);

    if o.RerunPatternSearch
        rerun=sim_k5p2_ps_recheck(B,R,casePlan,taskPlan,F, ...
            prodPath,i13dPath,outDir,'ParallelWorkers',o.ParallelWorkers);
        K52.optionalRecheck=rerun;
        K52.optimizationRerun=true;
    end
    save(fullfile(outDir,'K5p2_results.mat'),'K52','-v7.3');
    fprintf('COMMIT K5.2 COMPLETE | outputs: %s\n',outDir);
    fprintf('Figure: %s\n',K52.figurePaths.pdf);
    fprintf('Scientific method and all frozen original files remain unchanged.\n');
end


function [prod,zhou,i13d,k5]=resolve_sources(o)
    prod=char(string(o.ProductionBundlePath));
    if isempty(prod)
        parent=fullfile(fso_result_utils.default_results_root(),'i11p1_production');
        prod=find_single(fullfile(parent,'I11p1_production_*','i11p1ProductionBundle.mat'));
    end
    if exist(prod,'file')~=2,error('K5p2:Production','Missing I11.1 production MAT.');end
    prodRoot=fileparts(prod);
    zhou=char(string(o.ZhouResultPath));
    if isempty(zhou)
        zhou=fullfile(prodRoot,'i13_reviewer12','zhou_baseline','i13b_zhou_baseline.mat');
    end
    i13d=char(string(o.I13DResultPath));
    if isempty(i13d)
        i13d=find_single(fullfile(prodRoot,'i13_reviewer12', ...
            'final_optimizer_benchmark','I13D_*','i13d_optimizer_benchmark.mat'));
    end
    k5=char(string(o.K5ResultPath));
    if isempty(k5)
        k5=find_single(fullfile(fso_result_utils.revision2_results_root(), ...
            'commit_k5_fig7','K5_*','Fig7_four_method_results.mat'));
    end
    for path={zhou,i13d,k5}
        if exist(path{1},'file')~=2
            error('K5p2:SourceFile','Missing verified source: %s',path{1});
        end
    end
end

function p=find_single(pattern)
    a=dir(pattern);
    if numel(a)~=1
        error('K5p2:SourceDiscovery', ...
            'Found %d candidate files matching %s. Supply exact path explicitly.', ...
            numel(a),pattern);
    end
    p=fullfile(a(1).folder,a(1).name);
end

function verify_frozen(F,C,D)
    assert(strcmp(C.freezeVersion,'I11.1-v1') && strcmp(F.version,'I13D-v1') && ...
        strcmp(D.version,'I13D-v1') && ...
        isequal(F.MVec,[4 8 16 32]) && ...
        isequal(F.SNRVec,[20 30]) && ...
        isequal(F.TurbulenceVec,[0 0.3]) && ...
        isequal(F.Replicates,[1 2]) && ...
        isequal(F.RestartsByM,[6 6 8 8]) && ...
        isequal(F.MaxIterByM,[4000 4000 4000 8000]) && ...
        isequal(F.MaxEvalsPerStart,[4001 4001 4001 8001]), ...
        'K5p2:FrozenGrid','I13D original grid/budgets changed.');
    assert(C.P_avg==1 && F.P_avg==1 && D.P_avg==1 && ...
        F.minGap==0.01 && D.minGap==0.01 && C.minGap==0.01 && ...
        F.pinZero && D.pinZero && C.pinZero && ...
        F.tieToleranceBits==D.tieToleranceBits && ...
        isequal(F.PatternSearch,D.PatternSearch) && ...
        strcmp(C.FastFadingMethod,'logh') && ...
        strcmp(C.YGridMode,'candidate-adaptive') && ...
        C.YGridScale==1.1 && C.loghSpanSigma==8 && ...
        C.loghYBlockSize==512 && C.historyEvery==50, ...
        'K5p2:FrozenPolicy','Scientific constraints/evaluator changed.');
end

function outDir=prepare_directory(o)
    root=fso_result_utils.revision2_results_root();
    if strlength(string(o.ResumeDirectory))>0
        outDir=char(string(o.ResumeDirectory));
        if exist(outDir,'dir')~=7,error('K5p2:Resume','Missing directory.');end
    elseif strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
        if exist(outDir,'dir')==7
            error('K5p2:OutputExists','New directory must not exist; use ResumeDirectory.');
        end
    else
        stamp=char(datetime('now','TimeZone','UTC', ...
            'Format','yyyyMMdd''T''HHmmssSSS''Z'''));
        outDir=fullfile(root,'commit_k5p2_fig7',['K5p2_' stamp]);
    end
    croot=canonical(root);cpath=canonical(outDir);
    if ~(strcmp(cpath,croot)||startsWith(cpath,[croot filesep]))
        error('K5p2:OutputRoot','Output must remain under results/ieee_revision2.');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('K5p2:MakeDir','%s',msg);end
    end
end

function tf=same_path(a,b),tf=strcmp(canonical(a),canonical(b));end
function s=canonical(a),s=char(java.io.File(char(string(a))).getCanonicalPath());end
function t=is_text(v),t=ischar(v)||(isstring(v)&&isscalar(v));end
function t=is_bool(v),t=islogical(v)&&isscalar(v);end
function t=is_positive_int(v),t=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1&&v==fix(v);end
