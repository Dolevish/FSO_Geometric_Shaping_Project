function X0 = i13_shared_start_points(cfg,xBaseline,masterSeed,nStarts)
%I13_SHARED_START_POINTS  Reproduce the SA restart initial points exactly.
%
% Used by I.13 optimizer comparison so SA and the alternative optimizer see
% the same raw starts.  The construction matches sa_multistart/sa_restart_task:
%   restart 1: uniform-PAM baseline
%   restart k>=2: baseline + 0.2*N(0,I), Threefry substream k.
%
% Returned columns are RAW starts. Each optimizer is responsible for using
% the same canonical projection/constraints before its first evaluation.

    validateattributes(nStarts,{'numeric'},{'scalar','integer','positive','finite'});
    masterSeed=uint32(masterSeed);
    xBaseline=xBaseline(:);
    if numel(xBaseline)~=cfg.M
        error('i13_shared_start_points:Cardinality','Baseline cardinality mismatch.');
    end

    X0=zeros(cfg.M,nStarts);
    X0(:,1)=xBaseline;
    if nStarts==1,return;end

    s=RandStream('Threefry','Seed',double(masterSeed));
    for k=2:nStarts
        s.Substream=k;
        X0(:,k)=xBaseline+0.2*randn(s,cfg.M,1);
    end
end
