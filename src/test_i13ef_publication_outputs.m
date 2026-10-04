function report = test_i13ef_publication_outputs()
%TEST_I13EF_PUBLICATION_OUTPUTS  Table/report regression without graphics.

    fprintf('\n============================================================\n');
    fprintf('I.13E-F publication-output regression\n');
    fprintf('============================================================\n');

    result=synthetic_result();
    outDir=tempname;
    mkdir(outDir);
    cleanup=onCleanup(@() cleanup_dir(outDir)); %#ok<NASGU>

    O=make_i13ef_reviewer12_publication_outputs(result, ...
        'OutputDirectory',outDir,'MakeFigures',false);

    assert(height(O.pointwise)==16);
    assert(height(O.optimizerSummaryByM)==4);
    assert(height(O.turbulentComparison)==8);
    assert(height(O.turbulentGainSummaryByM)==4);
    assert(exist(O.paths.pointwise,'file')==2);
    assert(exist(O.paths.optimizerByM,'file')==2);
    assert(exist(O.paths.turbulent,'file')==2);
    assert(exist(O.paths.turbulentGainByM,'file')==2);
    assert(exist(O.paths.summaryText,'file')==2);

    assert(O.summary.nPhysicalPoints==16);
    assert(O.summary.nOptimizationCases==32);
    assert(abs(O.summary.tieToleranceBits-1e-3)<eps);

    report=struct('passed',true,'nPointwise',height(O.pointwise), ...
        'nTurbulent',height(O.turbulentComparison));
    fprintf('PASS: publication tables/report generated from synthetic frozen results.\n');
end


function result=synthetic_result()
    MVec=[4 8 16 32];
    SNRVec=[20 30];
    sigVec=[0 0.3];

    n=16;
    M=zeros(n,1); SNRdB=zeros(n,1); SigmaR2=zeros(n,1);
    PAMValidated=zeros(n,1); ZhouValidated=zeros(n,1);
    SAMean=zeros(n,1); PSMean=zeros(n,1); SAStd=zeros(n,1); PSStd=zeros(n,1);
    DeltaRep1=zeros(n,1); DeltaRep2=zeros(n,1); SAMinusPSMean=zeros(n,1);
    Outcome=strings(n,1);

    q=0;
    for m=MVec
        for snr=SNRVec
            for sig=sigVec
                q=q+1;
                M(q)=m; SNRdB(q)=snr; SigmaR2(q)=sig;
                PAMValidated(q)=0.5+0.05*log2(m)+0.01*snr-0.2*sig;
                ZhouValidated(q)=PAMValidated(q)+0.05+0.1*sig;
                PSMean(q)=ZhouValidated(q)+0.03+0.05*sig;
                SAMean(q)=PSMean(q)+0.002*log2(m)+0.01*sig;
                SAStd(q)=1e-3;
                PSStd(q)=2e-3;
                DeltaRep1(q)=SAMean(q)-PSMean(q)-2e-4;
                DeltaRep2(q)=SAMean(q)-PSMean(q)+2e-4;
                SAMinusPSMean(q)=mean([DeltaRep1(q) DeltaRep2(q)]);
                if abs(SAMinusPSMean(q))<=1e-3
                    Outcome(q)="Equivalent";
                elseif SAMinusPSMean(q)>0
                    Outcome(q)="SA better";
                else
                    Outcome(q)="Pattern Search better";
                end
            end
        end
    end

    physicalSummary=table(M,SNRdB,SigmaR2,PAMValidated,ZhouValidated, ...
        SAMean,PSMean,SAStd,PSStd,DeltaRep1,DeltaRep2,SAMinusPSMean,Outcome);

    nc=32;
    M=zeros(nc,1); PSRestartStd=zeros(nc,1);
    PSTotalFunctionEvaluations=zeros(nc,1); PSMaxEvalOverrun=zeros(nc,1);
    q=0;
    for m=MVec
        for snr=SNRVec
            for sig=sigVec
                for rep=1:2
                    q=q+1;
                    M(q)=m;
                    PSRestartStd(q)=0.01*log2(m);
                    PSTotalFunctionEvaluations(q)=1000;
                    PSMaxEvalOverrun(q)=0;
                end
            end
        end
    end
    caseSummary=table(M,PSRestartStd,PSTotalFunctionEvaluations,PSMaxEvalOverrun);

    config=struct('tieToleranceBits',1e-3,'totalMaxFunctionEvaluations',40000);
    metrics=struct('totalPatternSearchFunctionEvaluations',sum(PSTotalFunctionEvaluations));

    result=struct('physicalSummary',physicalSummary,'caseSummary',caseSummary, ...
        'config',config,'metrics',metrics,'outputDirectory',tempdir);
end


function cleanup_dir(path)
    if exist(path,'dir')==7
        try,rmdir(path,'s');catch,end
    end
end
