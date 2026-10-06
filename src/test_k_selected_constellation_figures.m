function report = test_k_selected_constellation_figures()
%TEST_K_SELECTED_CONSTELLATION_FIGURES
% Synthetic regression for Commit-K four-panel constellation figures.
% No optimizer or AMI evaluator is called.

    fprintf('\n============================================================\n');
    fprintf('COMMIT K - selected constellation figure regression\n');
    fprintf('============================================================\n');

    gaps=[0 0.005 0.01 0.02 0.05];
    MVec=[16 32];
    rows=cell(numel(gaps)*numel(MVec)*2,1);
    q=0;

    for M=MVec
        for d=gaps
            for rep=1:2
                q=q+1;
                base=(0:M-1).';
                x=base.^1.12;
                x=x/mean(x);
                if d>0
                    % Synthetic, visibly distinct geometry while preserving
                    % a positive increasing vector.
                    x=x + d*(0:M-1).';
                    x=x/mean(x);
                end
                rows{q}=struct( ...
                    'M',M,'SNRdB',20,'SigmaR2',0.2,'MinGap',d, ...
                    'Replicate',rep,'WinnerValidatedAMI', ...
                    1+0.001*M-0.3*d+0.002*rep, ...
                    'WinnerConstellation',x);
            end
        end
    end

    R=struct2table(vertcat(rows{:}));
    synthetic=struct('replicateSummary',R);

    outDir=tempname;
    mkdir(outDir);
    cleanup=onCleanup(@() cleanup_dir(outDir)); %#ok<NASGU>

    F=plot_k_selected_constellations(synthetic, ...
        'OutputDirectory',outDir,'CloseFigures',true);

    assert(height(F.plotSummary)==8);
    assert(all(F.plotSummary.SelectedReplicate==2));

    for M=MVec
        key=sprintf('M%d',M);
        assert(isfield(F.paths,key));
        assert(exist(F.paths.(key).pdf,'file')==2);
        assert(exist(F.paths.(key).png,'file')==2);
        assert(exist(F.paths.(key).fig,'file')==2);
    end

    csvPath=fullfile(outDir,'K_selected_constellation_plot_summary.csv');
    assert(exist(csvPath,'file')==2);

    fprintf('PASS: two 4-panel figures and summary CSV generated.\n');
    report=struct('passed',true,'nRows',height(F.plotSummary));
end


function cleanup_dir(path)
    if exist(path,'dir')==7
        try,rmdir(path,'s');catch,end
    end
end
