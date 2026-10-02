function report = test_i12_production_analysis()
%TEST_I12_PRODUCTION_ANALYSIS  Structural regression for I.12A analysis.
%
% This test builds a synthetic full-size I.11.1 bundle, writes it to a
% temporary directory, runs the analyzer, and checks aggregation semantics.
% It does not test the numerical AMI evaluators.

    fprintf('\n============================================================\n');
    fprintf('Commit I.12A - post-production analysis regression suite\n');
    fprintf('============================================================\n');

    tests={@test_full_synthetic_bundle};
    names={'Full 192-case synthetic bundle aggregates to 96 physical points'};

    passed=false(numel(tests),1); elapsed=nan(numel(tests),1); messages=strings(numel(tests),1);
    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k}); t=tic;
        try
            tests{k}(); passed(k)=true; messages(k)="PASS";
            elapsed(k)=toc(t); fprintf('PASS (%.2fs)\n',elapsed(k));
        catch ME
            elapsed(k)=toc(t); messages(k)=string(ME.message);
            fprintf('FAIL (%.2fs)\n  %s\n',elapsed(k),ME.message);
        end
    end

    report=struct('names',{names},'passed',passed,'elapsedSeconds',elapsed, ...
        'messages',messages,'nPassed',sum(passed),'nTotal',numel(tests));
    fprintf('\nPassed %d/%d tests in %.2fs.\n',report.nPassed,report.nTotal,sum(elapsed,'omitnan'));
    if ~all(passed)
        error('test_i12_production_analysis:Failure','I.12A analysis regression suite failed.');
    end
end

function test_full_synthetic_bundle()
    F=i11p1_production_config();
    root=tempname; mkdir(root);
    cleanup=onCleanup(@() cleanup_dir(root)); %#ok<NASGU>
    path=fullfile(root,'i11p1ProductionBundle.mat');

    [bundle,expected]=synthetic_bundle(F,root);
    save(path,'bundle','-v7.3');

    A=analyze_i11p1_production_results( ...
        'ProductionBundlePath',path,'OutputDirectory',fullfile(root,'analysis'),'WriteCSV',false);

    assert(height(A.caseTable)==192);
    assert(height(A.physicalTable)==96);
    assert(height(A.tableI_M8_SNR20)==4);
    assert(A.metrics.nSelectionChanges==expected.nSelectionChanges);
    assert(abs(A.metrics.maxFastValidatorGap-expected.maxGap)<1e-12);

    p=A.physicalTable(A.physicalTable.M==8 & ...
        abs(A.physicalTable.SigmaX2-0.2)<1e-12 & ...
        abs(A.physicalTable.SNRdB-20)<1e-12,:);
    assert(height(p)==1);
    assert(abs(p.GSMeanValidated-expected.specialMean)<1e-12);
    assert(abs(p.ReplicateRangeBits-expected.specialRange)<1e-12);
    assert(p.RepresentativeReplicate==1); % exact mean tie -> lower replicate

    assert(height(A.fig4_M8_SNR20)==4);
    assert(all(A.fig4_M8_SNR20.RepresentativeReplicate==1));
end

function [bundle,expected]=synthetic_bundle(F,root)
    rows=cell(F.nOptimizationCases,1); q=0;
    selectionCount=0; maxGap=0;
    specialMean=NaN; specialRange=NaN;

    for iM=1:numel(F.MVec)
        M=F.MVec(iM);
        [nStarts,~]=i11_restart_policy(M);
        dt=0.01;
        if (M==16 && any(abs(F.TurbulenceVec-[0.2 0.3]')<1e-12)) %#ok<NBRAK>
            % The exact synthetic dt values are not scientifically relevant.
        end
        for iSig=1:numel(F.TurbulenceVec)
            sig=F.TurbulenceVec(iSig);
            for iSNR=1:numel(F.SNRVec)
                snr=F.SNRVec(iSNR);
                pam=0.2+0.01*M+0.02*snr-0.1*sig;
                for rep=1:F.nReplicates
                    q=q+1;
                    gain=0.05+0.001*M+0.002*snr+0.01*sig+0.004*(rep-1);
                    ami=pam+gain;
                    gap=1e-5*(1+mod(q,7));
                    changed=mod(q,17)==0;
                    selectionCount=selectionCount+double(changed);
                    maxGap=max(maxGap,gap);
                    rows{q}=struct('M',M,'SigmaX2',sig,'SNRdB',snr, ...
                        'Replicate',rep,'NStarts',nStarts,'Seed',q, ...
                        'SelectedDt',dt,'AMIValidated',ami,'AMIFast',ami+gap/2, ...
                        'GainBits',gain,'K99',1000+q,'RestartValidatedStd',0.01+1e-5*q, ...
                        'MaxAbsFastValidationGap',gap,'SelectionChanged',changed, ...
                        'RestartTaskMaxSeconds',1,'CheckpointReusedRestarts',0);
                end
                if M==8 && abs(sig-0.2)<1e-12 && snr==20
                    g1=pam+(0.05+0.001*M+0.002*snr+0.01*sig);
                    g2=g1+0.004;
                    specialMean=(g1+g2)/2; specialRange=g2-g1;
                end
            end
        end
    end
    T=struct2table(vertcat(rows{:}));

    standard=make_stage([4 8 16],F,T);
    m32=make_stage(32,F,T);
    bundle=struct();
    bundle.productionFreeze=F;
    bundle.stagePlan=i11p1_stage_plan(F);
    bundle.run=struct('outputDirectory',root,'resumeRequested',false,'manifestPath','');
    bundle.stages=struct('m32',m32,'standard',standard);
    bundle.summary=struct('caseTable',T);

    expected=struct('nSelectionChanges',selectionCount,'maxGap',maxGap, ...
        'specialMean',specialMean,'specialRange',specialRange);
end

function stage=make_stage(MVec,F,T)
    sig=F.TurbulenceVec; snr=F.SNRVec; reps=1:F.nReplicates;
    nM=numel(MVec); nSig=numel(sig); nSNR=numel(snr); nRep=numel(reps);
    xWinner=cell(nM,nSig,nSNR,nRep); baseline=cell(nM,nSig,nSNR);

    for iM=1:nM
        M=MVec(iM);
        x=linspace(0,2,M).';
        for iSig=1:nSig
            for iSNR=1:nSNR
                baseline{iM,iSig,iSNR}=struct('x',x,'amiValidated',0);
                for r=1:nRep
                    % Small deterministic geometry perturbation with preserved mean.
                    d=(r-1)*1e-3*linspace(-1,1,M).';
                    xWinner{iM,iSig,iSNR,r}=x+d-mean(d);
                end
            end
        end
    end

    stage=struct();
    stage.axes=struct('M',MVec,'P_avg',1,'SNR_dB',snr,'sigma_X_sq',sig,'replicate',reps);
    stage.baseline=baseline;
    stage.summary=struct('xWinner',{xWinner});
    stage.parallel=struct('schedulerWallSeconds',0);
    stage.checkpoint=struct('nReused',0,'nExecuted',height(T));
end

function cleanup_dir(path)
    if exist(path,'dir')==7
        try,rmdir(path,'s');catch,end
    end
end
