classdef fso_result_utils
%FSO_RESULT_UTILS  Reproducible result/metadata helpers for IEEE revision.

methods(Static)

function root = default_results_root()
%DEFAULT_RESULTS_ROOT  Repository-local results directory.
    srcDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(srcDir);
    root = fullfile(repoRoot, 'results', 'ieee_revision');
end


function meta = run_metadata(sourceScript)
%RUN_METADATA  Common metadata stored with every experiment bundle.
    if nargin < 1 || isempty(sourceScript), sourceScript = 'unknown'; end

    meta = struct();
    meta.schemaVersion = 'FSO-GS-IEEE-revision-v1';
    meta.sourceScript  = char(sourceScript);
    meta.createdUTC    = char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
    meta.matlabVersion = version;
    meta.gitCommit     = 'unknown';

    srcDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(srcDir);
    try
        [status, txt] = system(sprintf('git -C "%s" rev-parse HEAD', repoRoot));
        if status == 0
            txt = strtrim(txt);
            if ~isempty(txt), meta.gitCommit = txt; end
        end
    catch
        % Metadata collection must never prevent a simulation from running.
    end
end


function snap = config_snapshot(cfg)
%CONFIG_SNAPSHOT  Serializable configuration without function handles/grid.
    snap = cfg;

    if isfield(snap, 'y_grid')
        yg = snap.y_grid;
        snap.yGridInfo = struct( ...
            'numPoints', numel(yg), ...
            'yMin', yg(1), ...
            'yMax', yg(end));
        snap = rmfield(snap, 'y_grid');
    end

    fieldsToRemove = {'AMI_Evaluator', 'AMI_Validator'};
    for k = 1:numel(fieldsToRemove)
        if isfield(snap, fieldsToRemove{k})
            snap = rmfield(snap, fieldsToRemove{k});
        end
    end
end


function seed = case_seed(baseSeed, M, SNR_dB, sigma_X_sq)
%CASE_SEED  Deterministic uint32-compatible seed for one channel case.
% The mapping is stable and independent of loop/worker scheduling.
    validateattributes(baseSeed, {'numeric'}, {'scalar','real','finite','nonnegative'});

    s = uint64(floor(double(baseSeed)));
    s = s + uint64(M) * uint64(1000003);
    s = s + uint64(round((double(SNR_dB) + 200) * 100)) * uint64(9176);
    s = s + uint64(round(double(sigma_X_sq) * 1e6)) * uint64(131);

    modulus = uint64(4294967295); % 2^32-1
    seed = double(mod(s, modulus));
end


function xmax = feasible_xmax_bound(M, P_avg, minGap, pinZero)
%FEASIBLE_XMAX_BOUND  Conservative maximum level under mean-intensity power.
% For pinZero=true and fixed d_min, the exact largest feasible final level is
% obtained when all free power is placed in the last residual coordinate.
    if nargin < 4 || isempty(pinZero), pinZero = true; end

    validateattributes(M, {'numeric'}, {'scalar','integer','>=',2});
    validateattributes(P_avg, {'numeric'}, {'scalar','real','positive','finite'});
    validateattributes(minGap, {'numeric'}, {'scalar','real','nonnegative','finite'});

    if pinZero
        xmax = M * P_avg - 0.5 * minGap * (M-1) * (M-2);
    else
        xmax = M * P_avg;
    end

    xmax = max(xmax, 2 * P_avg);
end


function name = case_filename(M, SNR_dB, sigma_X_sq)
%CASE_FILENAME  Stable filesystem-safe case filename.
    sigToken = strrep(sprintf('%.4f', sigma_X_sq), '.', 'p');
    snrToken = strrep(sprintf('%+.2f', SNR_dB), '.', 'p');
    snrToken = strrep(snrToken, '+', 'p');
    snrToken = strrep(snrToken, '-', 'm');
    name = sprintf('case_M%d_SNR%s_sig%s.mat', M, snrToken, sigToken);
end

end % methods(Static)
end % classdef
