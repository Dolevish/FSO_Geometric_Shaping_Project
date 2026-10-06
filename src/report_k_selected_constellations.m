function report = report_k_selected_constellations(resultOrPath,varargin)
%REPORT_K_SELECTED_CONSTELLATIONS  Print Commit-K selected-geometry tables.
%
% For every physical Commit-K point (M,d_min), this post-processor selects
% the replicate with the LARGER validated AMI. It does not rerun SA or the
% AMI evaluators.
%
% Usage:
%   report = report_k_selected_constellations(Rk);
%
% or, after restarting MATLAB:
%   report = report_k_selected_constellations( ...
%       '/path/to/K_min_gap_results.mat');
%
% Printed outputs:
%   1) Selected constellation levels for all 10 (M,d_min) points.
%   2) Compact table: M, d_min, d_actual, validated AMI.
%   3) Detailed diagnostic table including gap slack, binding status,
%      number of active adjacent pairs, and geometry statistics.
%
% Binding definition:
%   abs(d_actual-d_min) <= BindingTolerance
% and NumActivePairs counts adjacent gaps satisfying
%   abs((x_{i+1}-x_i)-d_min) <= BindingTolerance.
%
% No global-optimum claim is implied.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'BindingTolerance',1e-4,@nonnegative_scalar);
    addParameter(p,'DisplayPrecision',8,@positive_integer);
    parse(p,varargin{:}); o=p.Results;

    result=load_result(resultOrPath);

    if ~isfield(result,'replicateSummary') || ~istable(result.replicateSummary)
        error('report_k_selected_constellations:MissingReplicateSummary', ...
            'Commit-K result must contain replicateSummary.');
    end

    R=result.replicateSummary;
    required={'M','SNRdB','SigmaR2','MinGap','Replicate','WinnerRestart', ...
        'WinnerFastAMI','WinnerValidatedAMI','WinnerConstellation'};
    missing=setdiff(required,R.Properties.VariableNames);
    if ~isempty(missing)
        error('report_k_selected_constellations:MissingFields', ...
            'replicateSummary is missing: %s',strjoin(missing,', '));
    end

    keys=unique(R(:,{'M','SNRdB','SigmaR2','MinGap'}),'rows','sorted');
    n=height(keys);

    M=zeros(n,1);
    SNRdB=zeros(n,1);
    SigmaR2=zeros(n,1);
    d_min=zeros(n,1);
    SelectedReplicate=zeros(n,1);
    WinnerRestart=zeros(n,1);
    FastAMI=zeros(n,1);
    ValidatedAMI=zeros(n,1);
    d_actual=zeros(n,1);
    GapSlack=zeros(n,1);
    BindingFlag=false(n,1);
    NumActivePairs=zeros(n,1);
    MaxLevel=zeros(n,1);
    MeanAdjacentGap=zeros(n,1);
    MedianAdjacentGap=zeros(n,1);
    ConstellationLevels=strings(n,1);
    LevelDiffs=strings(n,1);
    LevelVector=cell(n,1);
    DiffVector=cell(n,1);

    for k=1:n
        idx=R.M==keys.M(k) & ...
            abs(R.SNRdB-keys.SNRdB(k))<1e-12 & ...
            abs(R.SigmaR2-keys.SigmaR2(k))<1e-12 & ...
            abs(R.MinGap-keys.MinGap(k))<1e-12;

        G=R(idx,:);
        if isempty(G)
            error('report_k_selected_constellations:MissingGroup', ...
                'No replicate rows found for M=%d, d_min=%.6g.', ...
                keys.M(k),keys.MinGap(k));
        end

        % Deterministic tie-break: largest validated AMI, then lowest
        % replicate number.
        G=sortrows(G,{'WinnerValidatedAMI','Replicate'},{'descend','ascend'});
        W=G(1,:);

        x=sort(double(W.WinnerConstellation{1}(:)));
        dx=diff(x);
        if isempty(dx)
            actual=NaN;
            meanGap=NaN;
            medianGap=NaN;
            nActive=0;
            binding=false;
        else
            actual=min(dx);
            meanGap=mean(dx);
            medianGap=median(dx);
            nActive=sum(abs(dx-double(W.MinGap))<=double(o.BindingTolerance));
            binding=abs(actual-double(W.MinGap))<=double(o.BindingTolerance);
        end

        M(k)=W.M;
        SNRdB(k)=W.SNRdB;
        SigmaR2(k)=W.SigmaR2;
        d_min(k)=W.MinGap;
        SelectedReplicate(k)=W.Replicate;
        WinnerRestart(k)=W.WinnerRestart;
        FastAMI(k)=W.WinnerFastAMI;
        ValidatedAMI(k)=W.WinnerValidatedAMI;
        d_actual(k)=actual;
        GapSlack(k)=actual-W.MinGap;
        BindingFlag(k)=binding;
        NumActivePairs(k)=nActive;
        MaxLevel(k)=max(x);
        MeanAdjacentGap(k)=meanGap;
        MedianAdjacentGap(k)=medianGap;
        ConstellationLevels(k)=string(mat2str(x.',double(o.DisplayPrecision)));
        LevelDiffs(k)=string(mat2str(dx.',double(o.DisplayPrecision)));
        LevelVector{k}=x(:).';
        DiffVector{k}=dx(:).';
    end

    % Table 1: complete selected constellation vectors.
    constellationTable=table(M,d_min,SelectedReplicate,ValidatedAMI,ConstellationLevels);

    % Table 2: compact reviewer-facing spacing comparison requested by user.
    compactTable=table(M,d_min,d_actual,ValidatedAMI);

    % Table 3: detailed selected-replicate diagnostics.
    detailedTable=table(M,SNRdB,SigmaR2,d_min,SelectedReplicate,WinnerRestart, ...
        ValidatedAMI,FastAMI,d_actual,GapSlack,BindingFlag,NumActivePairs, ...
        MaxLevel,MeanAdjacentGap,MedianAdjacentGap,ConstellationLevels,LevelDiffs);

    fprintf('\n============================================================\n');
    fprintf('COMMIT K - SELECTED CONSTELLATIONS\n');
    fprintf('Selection rule: highest validated AMI replicate for each (M,d_min).\n');
    fprintf('Binding tolerance: %.3e\n',o.BindingTolerance);
    fprintf('No optimization or AMI evaluation was rerun.\n');
    fprintf('============================================================\n');

    fprintf('\nTABLE 1 - FINAL CONSTELLATION LEVELS\n');
    disp(constellationTable);

    fprintf('\nFull constellation vectors (untruncated):\n');
    for k=1:n
        fprintf('M=%d | d_min=%.3f | rep=%d | validated AMI=%.8f | x=%s\n', ...
            M(k),d_min(k),SelectedReplicate(k),ValidatedAMI(k), ...
            char(ConstellationLevels(k)));
    end

    fprintf('\nTABLE 2 - IMPOSED VS ACTUAL MINIMUM SPACING\n');
    disp(compactTable);

    fprintf('\nTABLE 3 - DETAILED SELECTED-CONSTELLATION DIAGNOSTICS\n');
    disp(detailedTable);

    fprintf('\nDetailed level differences (untruncated):\n');
    for k=1:n
        fprintf(['M=%d | d_min=%.3f | d_actual=%.8f | slack=%+.8f | ' ...
            'binding=%d | active_pairs=%d | dx=%s\n'], ...
            M(k),d_min(k),d_actual(k),GapSlack(k),BindingFlag(k), ...
            NumActivePairs(k),char(LevelDiffs(k)));
    end
    fprintf('============================================================\n');

    report=struct();
    report.selectionRule='highest validated AMI replicate per (M,d_min)';
    report.bindingTolerance=double(o.BindingTolerance);
    report.constellationTable=constellationTable;
    report.compactTable=compactTable;
    report.detailedTable=detailedTable;
    report.levelVectors=LevelVector;
    report.diffVectors=DiffVector;
    report.globalOptimumClaim=false;
end


function result=load_result(v)
    if isstruct(v)
        result=v;
        return;
    end

    if ~(ischar(v)||(isstring(v)&&isscalar(v)))
        error('report_k_selected_constellations:Input', ...
            'Input must be a Commit-K result struct or MAT-file path.');
    end

    path=char(string(v));
    if exist(path,'file')~=2
        error('report_k_selected_constellations:MissingFile', ...
            'MAT file does not exist: %s',path);
    end

    S=load(path);
    if isfield(S,'result')
        result=S.result;
    elseif isfield(S,'Rk')
        result=S.Rk;
    else
        error('report_k_selected_constellations:MissingResult', ...
            'MAT file must contain variable result or Rk.');
    end
end


function tf=nonnegative_scalar(v)
    tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;
end

function tf=positive_integer(v)
    tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0&&mod(v,1)==0;
end
