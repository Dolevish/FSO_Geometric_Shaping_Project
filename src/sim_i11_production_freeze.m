function bundle = sim_i11_production_freeze(varargin)
%SIM_I11_PRODUCTION_FREEZE  Compatibility alias for the current production freeze.
%
% I.11-v1 used one fixed 4000-iteration budget. The canonical production
% profile was superseded by I.11.1, which keeps the frozen restart map but
% uses 8000 iterations/restart for M=32 and 4000 for M=4,8,16.
%
% This alias intentionally forwards to sim_i11p1_production_freeze so old
% launch commands cannot accidentally start the superseded production run.

    bundle=sim_i11p1_production_freeze(varargin{:});
end
