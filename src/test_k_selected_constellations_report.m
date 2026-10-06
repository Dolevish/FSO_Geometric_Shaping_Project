function report = test_k_selected_constellations_report()
%TEST_K_SELECTED_CONSTELLATIONS_REPORT  Synthetic regression for K report.
%
% No optimization/evaluator calls are made.

    fprintf('\n============================================================\n');
    fprintf('COMMIT K - selected-constellation report regression\n');
    fprintf('============================================================\n');

    rows=cell(4,1);
    q=0;
    for d=[0 0.01]
        for rep=1:2
            q=q+1;
            x=[0;0.25+d;0.75+d;3-2*d];
            rows{q}=struct( ...
                'M',4,'SNRdB',20,'SigmaR2',0.2,'MinGap',d, ...
                'Replicate',rep,'WinnerRestart',rep, ...
                'WinnerFastAMI',1+0.01*rep-0.1*d, ...
                'WinnerValidatedAMI',1+0.02*rep-0.1*d, ...
                'WinnerConstellation',x);
        end
    end
    R=struct2table(vertcat(rows{:}));
    synthetic=struct('replicateSummary',R);

    out=evalc('Q=report_k_selected_constellations(synthetic);');

    assert(height(Q.constellationTable)==2);
    assert(height(Q.compactTable)==2);
    assert(height(Q.detailedTable)==2);
    assert(all(Q.constellationTable.SelectedReplicate==2));
    assert(all(abs(Q.compactTable.ValidatedAMI-[1.04;1.039])<1e-12));
    assert(contains(out,'TABLE 1 - FINAL CONSTELLATION LEVELS'));
    assert(contains(out,'TABLE 2 - IMPOSED VS ACTUAL MINIMUM SPACING'));
    assert(contains(out,'TABLE 3 - DETAILED SELECTED-CONSTELLATION DIAGNOSTICS'));

    fprintf('PASS: best replicate selected and all three requested tables printed.\n');

    report=struct('passed',true,'capturedOutput',out);
end
