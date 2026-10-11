function P=select_i11p1_best_validated_curves(A)
%SELECT_I11P1_BEST_VALIDATED_CURVES  Observed best validated AMI per channel point.
%
% The frozen I.11.1 case table has one row per (M, sigma_R^2, SNR, replicate).
% Its AMIValidated is the MAX validated score across ALL SA restarts for that
% replicate (see sa_finalize_restart_results). Thus max over replicates is
% exactly max over all saved validated production restarts for that point.
% Never use the FAST optimizer score for this selection.
%
% This function is read-only: no SA, AMI calls, channel sampling, or data edits.
% Existing GSMeanValidated and other historical summary columns are retained.

    if ~isstruct(A) || ~all(isfield(A,{'physicalTable','caseTable','productionFreeze'}))
        error('select_i11p1_best_validated_curves:BadAnalysis','Expected I.12A analysis.');
    end
    P=A.physicalTable;
    T=A.caseTable;
    requiredP={'M','SigmaX2','SNRdB','PAMValidated', ...
        'GSRep1Validated','GSRep2Validated','GSMeanValidated'};
    requiredT={'M','SigmaX2','SNRdB','Replicate','NStarts', ...
        'AMIValidated','GainBits'};
    if ~istable(P) || ~istable(T) || ...
            ~all(ismember(requiredP,P.Properties.VariableNames)) || ...
            ~all(ismember(requiredT,T.Properties.VariableNames))
        error('select_i11p1_best_validated_curves:MissingColumns', ...
            'Analysis must contain I.11.1 physical and replicate case tables.');
    end

    nRep=A.productionFreeze.nReplicates;
    if nRep~=2 || height(T)~=height(P)*nRep
        error('select_i11p1_best_validated_curves:CaseCount', ...
            'Expected exactly two saved replicates per physical point.');
    end

    n=height(P);
    bestAMI=nan(n,1);
    bestRep=nan(n,1);
    considered=nan(n,1);
    gainBest=nan(n,1);
    for k=1:n
        idx=(T.M==P.M(k)) & ...
            abs(T.SigmaX2-P.SigmaX2(k))<1e-12 & ...
            abs(T.SNRdB-P.SNRdB(k))<1e-12;
        R=sortrows(T(idx,:),'Replicate');
        if height(R)~=nRep || ~isequal(R.Replicate(:).',1:nRep)
            error('select_i11p1_best_validated_curves:MissingReplicate', ...
                'Missing or duplicated replicate at physical row %d.',k);
        end
        vals=R.AMIValidated(:);
        baselines=vals-R.GainBits(:);
        if any(~isfinite(vals)) || any(~isfinite(baselines)) || ...
                any(~isfinite(R.NStarts)) || any(R.NStarts<1)
            error('select_i11p1_best_validated_curves:NonFinite', ...
                'Invalid saved validated scores or restart counts, row %d.',k);
        end
        if abs(vals(1)-P.GSRep1Validated(k))>1e-9 || ...
                abs(vals(2)-P.GSRep2Validated(k))>1e-9 || ...
                max(abs(baselines-P.PAMValidated(k)))>1e-9 || ...
                abs(mean(vals)-P.GSMeanValidated(k))>1e-9
            error('select_i11p1_best_validated_curves:InconsistentSummary', ...
                'Frozen replicate and physical tables disagree, row %d.',k);
        end

        % sortrows puts replicate 1 first. max uses the first in a tie.
        [bestAMI(k),winner]=max(vals);
        bestRep(k)=R.Replicate(winner);
        considered(k)=sum(R.NStarts);
        gainBest(k)=bestAMI(k)-P.PAMValidated(k);
    end

    P.GSBestValidated=bestAMI;
    P.GSBestReplicate=bestRep;
    P.GSCurveValidated=bestAMI;
    P.GainBestBits=gainBest;
    P.NRestartsConsidered=considered;
end
