function out=make_ieee_revision_best_validated_legacy_curves(varargin)
%MAKE_IEEE_REVISION_BEST_VALIDATED_LEGACY_CURVES
% Replot the two I.12C paper performance figures using the observed highest
% VALIDATED AMI of all I.11.1 production SA restarts at each channel point.
% The M=8 and all-M graph styles/axes/legends are inherited unchanged from
% make_ieee_revision_figures_legacy_style. No SA/AMI evaluation is performed.
%
%   out=make_ieee_revision_best_validated_legacy_curves();
%   out=make_ieee_revision_best_validated_legacy_curves( ...
%       'ProductionBundlePath','/path/to/i11p1ProductionBundle.mat', ...
%       'CloseFigures',true);
%
% Writes PDF/PNG/FIG plus selected-point CSV and provenance manifest to a
% NEW directory. The original mean-curve figures and Fig. 4 stay untouched.
    out=make_ieee_revision_figures_legacy_style( ...
        'PerformanceAggregation','best-validated',varargin{:});
end
