function [out,results] = sa_finalize_restart_results(results,taskLabel)
%SA_FINALIZE_RESTART_RESULTS  Rank independently executed restart results.
%
% This is the deterministic reduction stage used by the I.10.5A flattened
% restart scheduler.  It reproduces the validated-selection logic of
% sa_multistart after every restart has already completed and been validated.

    if nargin < 2 || isempty(taskLabel), taskLabel='[SA]'; end
    if isstring(taskLabel), taskLabel=char(taskLabel); end
    if isempty(results)
        error('sa_finalize_restart_results:EmptyResults','At least one restart result is required.');
    end

    [~,ord] = sort([results.startIndex]);
    results = results(ord);
    nStarts = numel(results);

    expected = 1:nStarts;
    if ~isequal([results.startIndex],expected)
        error('sa_finalize_restart_results:NonContiguousStarts', ...
            'Restart indices must be exactly 1:nStarts for final selection.');
    end

    allFast = [results.bestMIFast];
    allVal = [results.bestMIValidated];
    [~,fastOrder] = sort(allFast,'descend');
    [~,valOrder] = sort(allVal,'descend');
    for r=1:nStarts
        results(fastOrder(r)).fastRank=r;
        results(valOrder(r)).validatedRank=r;
    end

    [bestMIFast,idxFast]=max(allFast);
    [bestMIVal,idxVal]=max(allVal);
    selectionChanged=idxFast~=idxVal;
    validatedAdvantage=bestMIVal-results(idxFast).bestMIValidated;
    absGaps=abs([results.validationGap]);

    out=struct();
    out.masterSeed=results(1).masterSeed;
    out.nStarts=nStarts;

    out.bestStartFast=idxFast;
    out.bestMIFast=bestMIFast;
    out.bestXFast=results(idxFast).x_best(:);
    out.fastWinnerValidatedMI=results(idxFast).bestMIValidated;

    out.bestStartValidated=idxVal;
    out.bestMIValidated=bestMIVal;
    out.bestXValidated=results(idxVal).x_best(:);
    out.validatedWinnerFastMI=results(idxVal).bestMIFast;

    out.selectionChanged=selectionChanged;
    out.validatedAdvantageOverFastWinner=validatedAdvantage;
    out.meanAbsFastValidationGap=mean(absGaps,'omitnan');
    out.maxAbsFastValidationGap=max(absGaps);
    out.validationTotalRuntime=sum([results.validationRuntime]);
    out.totalSARuntime=sum([results.runtime]);
    out.totalRuntime=out.totalSARuntime+out.validationTotalRuntime;
    out.meanAcceptanceRate=mean([results.acceptanceRate],'omitnan');

    out.bestStart=out.bestStartValidated;
    out.bestMI=out.bestMIValidated;
    out.bestX=out.bestXValidated;
    out.bestRun=results(idxVal);

    fprintf('%s selection summary | fast winner=%d | validated winner=%d | changed=%d\n', ...
        taskLabel,idxFast,idxVal,selectionChanged);
    fprintf('%s final validated AMI=%.6f | fast-winner validated AMI=%.6f | improvement=%.3e\n', ...
        taskLabel,bestMIVal,results(idxFast).bestMIValidated,validatedAdvantage);
end
