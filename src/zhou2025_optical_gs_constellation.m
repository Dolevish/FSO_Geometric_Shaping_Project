function out = zhou2025_optical_gs_constellation(M,Pavg,varargin)
%ZHOU2025_OPTICAL_GS_CONSTELLATION  Literature GS baseline from Zhou et al.
%
%   out = zhou2025_optical_gs_constellation(M,Pavg)
%   out = zhou2025_optical_gs_constellation(...,'ExtraBits',2)
%
% Reproduces the optical-intensity geometric-shaping construction in:
%   S. Zhou et al., "Geometric Constellation Shaping for Wireless Optical
%   Intensity Channels: An Information-Theoretic Approach,"
%   IEEE Communications Letters, vol. 29, no. 1, 2025.
%
% The implementation follows the paper's equations (5)-(19):
%   1) exponential-input quantiles and interval centroids c_m;
%   2) shifting/scaling to X_l with first level zero and unchanged mean;
%   3) finite-resolution regularization with b+n bits.
%
% This function is intentionally CHANNEL INDEPENDENT.  It does not use the
% FSO turbulence law, SNR, or our AMI optimizer.  That is precisely why it
% is useful as a literature GS baseline in Reviewer-1 analysis.
%
% Inputs
%   M       constellation size; power of two (project uses 4,8,16,32)
%   Pavg    mean optical intensity E > 0
%
% Name-value
%   ExtraBits  n in the paper's b+n-bit regularization (default 2)
%   Regularize true -> return the paper's regularized X_ell as out.x
%              false -> return shifted/scaled X_l as out.x
%
% Output fields
%   x                 selected baseline constellation
%   xCentroid         X_c
%   xShiftScale       X_l
%   xRegularized      X_ell
%   integerLevels     ell_m
%   deltaGrid         d from Eq. (17)
%   basicLevel        Delta from Eq. (19)
%   meanPower, minGap, PAPR
%   metadata          source/equation information
%
% No claim is made that this construction is optimal for the No-CSI FSO
% channel.  It is evaluated there only as an external literature baseline.

    p=inputParser; p.FunctionName=mfilename;
    addRequired(p,'M',@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=2&&mod(v,1)==0);
    addRequired(p,'Pavg',@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0);
    addParameter(p,'ExtraBits',2,@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0&&mod(v,1)==0);
    addParameter(p,'Regularize',true,@(v)islogical(v)&&isscalar(v));
    parse(p,M,Pavg,varargin{:}); o=p.Results;

    M=double(M); E=double(Pavg); n=double(o.ExtraBits);
    b=log2(M);
    if abs(b-round(b))>10*eps(max(1,b))
        error('zhou2025_optical_gs_constellation:PowerOfTwo', ...
            'The b+n-bit regularization requires M to be a power of two.');
    end
    b=round(b);

    % ------------------------------------------------------------------
    % Eqs. (5)-(8): equiprobable exponential intervals and centroids.
    % q_0=0, q_M=infinity, q_m=E ln(M/(M-m)).
    % ------------------------------------------------------------------
    q=nan(M+1,1);
    q(1)=0;
    for m=1:M-1
        q(m+1)=E*log(M/(M-m));
    end
    q(M+1)=Inf;

    c=zeros(M,1);
    for m=0:M-2
        qm=q(m+1);
        qmp1=q(m+2);
        c(m+1)=M*((E+qm)*exp(-qm/E) - (E+qmp1)*exp(-qmp1/E));
    end
    c(M)=E*(log(M)+1); % Eq. (8)

    % ------------------------------------------------------------------
    % Eqs. (13)-(14): shift c0 to zero, then scale to restore mean E.
    % ------------------------------------------------------------------
    scaleDen=(M-1)*log(M/(M-1));
    g=1/scaleDen;
    xShift=(c-c(1))*g;

    % Numerical consistency checks for the published construction.
    tol=1e-11*max(1,E);
    if abs(mean(c)-E)>tol
        error('zhou2025_optical_gs_constellation:CentroidMean', ...
            'Centroid construction mean mismatch: %.3e.',mean(c)-E);
    end
    if abs(mean(xShift)-E)>tol || abs(xShift(1))>tol
        error('zhou2025_optical_gs_constellation:ShiftScaleMean', ...
            'Shift/scale construction failed mean or zero-anchor check.');
    end

    % ------------------------------------------------------------------
    % Eqs. (17)-(19): b+n-bit regularization and mean-power fine tuning.
    % ------------------------------------------------------------------
    maxInteger=2^(b+n)-1;
    d=xShift(end)/maxInteger;
    ell=floor(xShift/d + 0.5);
    ell(1)=0; % exact zero anchor

    if mean(ell)<=0
        error('zhou2025_optical_gs_constellation:DegenerateRegularization', ...
            'Regularized integer levels have zero mean.');
    end
    Delta=E/mean(ell);
    xReg=double(ell)*Delta;

    if abs(mean(xReg)-E)>tol
        error('zhou2025_optical_gs_constellation:RegularizedMean', ...
            'Regularized construction mean mismatch.');
    end
    if any(diff(xReg)<-tol)
        error('zhou2025_optical_gs_constellation:Ordering', ...
            'Regularized levels are not nondecreasing.');
    end

    if o.Regularize
        x=xReg;
        variant='regularized-b-plus-n';
    else
        x=xShift;
        variant='shifted-scaled-continuous';
    end

    out=struct();
    out.M=M;
    out.Pavg=E;
    out.extraBits=n;
    out.bits=b+n;
    out.variant=variant;
    out.quantiles=q;
    out.xCentroid=c;
    out.shiftScaleFactor=g;
    out.xShiftScale=xShift;
    out.integerLevels=ell;
    out.deltaGrid=d;
    out.basicLevel=Delta;
    out.xRegularized=xReg;
    out.x=x(:);
    out.meanPower=mean(x);
    out.minGap=min(diff(x));
    out.PAPR=max(x)/mean(x);
    out.metadata=struct( ...
        'source','Zhou et al., IEEE Communications Letters, vol. 29, no. 1, 2025', ...
        'construction','Eqs. (5)-(19): centroid -> shift/scale -> regularization', ...
        'channelAdapted',false, ...
        'usesFSOTurbulence',false, ...
        'globalOptimumClaim',false);
end
