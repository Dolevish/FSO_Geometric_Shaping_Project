function report=test_best_validated_curve_selection()
%TEST_BEST_VALIDATED_CURVE_SELECTION  No MATLAB optimization or AMI calls.
% A permutation-invariant, synthetic-only unit test of max-validated choice.

    P=table([8;16],[0.1;0.3],[20;25],[1.2;1.8], ...
        [1.52;2.1],[1.55;2.1],[1.535;2.1], ...
        'VariableNames',{'M','SigmaX2','SNRdB','PAMValidated', ...
        'GSRep1Validated','GSRep2Validated','GSMeanValidated'});
    T=table([16;8;16;8],[0.3;0.1;0.3;0.1], ...
        [25;20;25;20],[2;2;1;1],[8;6;8;6], ...
        [2.1;1.55;2.1;1.52],[0.3;0.35;0.3;0.32], ...
        'VariableNames',{'M','SigmaX2','SNRdB','Replicate', ...
        'NStarts','AMIValidated','GainBits'});
    A=struct('physicalTable',P,'caseTable',T, ...
        'productionFreeze',struct('nReplicates',2));

    S=select_i11p1_best_validated_curves(A);
    assert(isequal(S.GSBestReplicate,[2;1]),'Best replicate/tie break incorrect.');
    assert(max(abs(S.GSBestValidated-[1.55;2.1]))<1e-12);
    assert(max(abs(S.GainBestBits-[0.35;0.3]))<1e-12);
    assert(isequal(S.NRestartsConsidered,[12;16]));
    assert(isequal(S.GSCurveValidated,S.GSBestValidated));
    assert(isequal(S.GSMeanValidated,P.GSMeanValidated), ...
        'The historical mean summary was modified.');

    A.caseTable=T([4 3 2 1],:);
    Z=select_i11p1_best_validated_curves(A);
    assert(isequal(S.GSBestValidated,Z.GSBestValidated) && ...
        isequal(S.GSBestReplicate,Z.GSBestReplicate), ...
        'Selection changed under row permutation.');

    A.caseTable=T(1:3,:);
    must_fail(@() select_i11p1_best_validated_curves(A));
    A.caseTable=T;
    A.physicalTable.GSRep1Validated(1)=999;
    must_fail(@() select_i11p1_best_validated_curves(A));

    report=struct('passed',true,'syntheticOnly',true,'nPhysicalPoints',height(P));
    fprintf('PASS: best-validated selection, ties, permutations, constraints on saved summaries.\n');
end

function must_fail(fn)
    try
        fn();
    catch
        return;
    end
    error('test_best_validated_curve_selection:ExpectedError', ...
        'Invalid source data unexpectedly passed validation.');
end
