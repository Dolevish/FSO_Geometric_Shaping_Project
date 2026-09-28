function report = test_core_invariants(varargin)
%TEST_CORE_INVARIANTS  Regression/invariant suite for the IEEE revision core.
%
%   report = test_core_invariants()
%   report = test_core_invariants('Name',Value,...)
%
% The suite is intentionally small enough to run before every research
% experiment, while covering the invariants introduced in Commits A-D:
%   * canonical SNR/configuration convention,
%   * exact mean-power + min-gap projection,
%   * projection idempotence,
%   * explicit pinZero semantics,
%   * infeasible min-gap rejection,
%   * Uniform-PAM baseline feasibility,
%   * fast/validated AMI consistency,
%   * deterministic case seeds,
%   * SA repeatability under a fixed seed,
%   * convergence-history accounting,
%   * final multi-start winner selected by validated AMI.
%
% Name-value options:
%   'ProjectionTrials'       random projection trials per M (default 20)
%   'ProjectionTolerance'    geometric/power tolerance (default 1e-10)
%   'AMITolerance'           allowed |fast-validation| gap (default 5e-3)
%   'ReproTolerance'         same-seed numerical tolerance (default 1e-10)
%   'saMaxIter'              short SA iterations for regression test (default 200)
%   'saNStarts'              SA restarts for regression test (default 2)
%   'seedInit'               deterministic SA seed (default 12345)
%   'Verbose'                print PASS/FAIL lines (default true)
%
% The function throws test_core_invariants:Failed if any test fails.
% It does NOT write simulation result files.

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'ProjectionTrials', 20, @(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p, 'ProjectionTolerance', 1e-10, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p, 'AMITolerance', 5e-3, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p, 'ReproTolerance', 1e-10, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p, 'saMaxIter', 200, @(v) isnumeric(v) && isscalar(v) && v>=1 && mod(v,1)==0);
    addParameter(p, 'saNStarts', 2, @(v) isnumeric(v) && isscalar(v) && v>=2 && mod(v,1)==0);
    addParameter(p, 'seedInit', 12345, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p, 'Verbose', true, @(v) islogical(v) && isscalar(v));
    parse(p, varargin{:});
    o = p.Results;

    testNames = {
        'Canonical configuration and SNR convention'
        'Uniform-PAM baseline feasibility'
        'Projection feasibility and idempotence'
        'Projection with pinZero=false'
        'Infeasible minimum-gap rejection'
        'Fast vs validated AMI consistency'
        'Deterministic case-seed mapping'
        'Same-seed SA repeatability and observability'
        'Validated multi-start winner consistency'
        };

    testFns = {
        @() test_config_and_snr(o)
        @() test_uniform_pam(o)
        @() test_projection(o)
        @() test_projection_unpinned(o)
        @() test_infeasible_gap()
        @() test_ami_consistency(o)
        @() test_case_seed()
        @() test_sa_repeatability(o)
        @() test_validated_winner(o)
        };

    nTests = numel(testNames);
    passed = false(nTests,1);
    elapsed = nan(nTests,1);
    messages = strings(nTests,1);

    if o.Verbose
        fprintf('\n============================================================\n');
        fprintf('FSO GS IEEE revision - core regression/invariant suite\n');
        fprintf('Projection trials/M=%d | SA=%d starts x %d iterations\n', ...
            o.ProjectionTrials, o.saNStarts, o.saMaxIter);
        fprintf('============================================================\n');
    end

    for k = 1:nTests
        if o.Verbose
            fprintf('[%02d/%02d] %s ... ', k, nTests, testNames{k});
        end
        t = tic;
        try
            testFns{k}();
            elapsed(k) = toc(t);
            passed(k) = true;
            messages(k) = "PASS";
            if o.Verbose
                fprintf('PASS (%.2fs)\n', elapsed(k));
            end
        catch ME
            elapsed(k) = toc(t);
            passed(k) = false;
            messages(k) = string(ME.message);
            if o.Verbose
                fprintf('FAIL (%.2fs)\n', elapsed(k));
                fprintf('        %s [%s]\n', ME.message, ME.identifier);
            end
        end
    end

    report = table(string(testNames), passed, elapsed, messages, ...
        'VariableNames', {'Test','Passed','ElapsedSeconds','Message'});

    if o.Verbose
        fprintf('\n------------------------- SUMMARY --------------------------\n');
        fprintf('Passed %d/%d tests in %.2fs.\n', nnz(passed), nTests, sum(elapsed,'omitnan'));
        if all(passed)
            fprintf('CORE INVARIANTS: PASS\n');
        else
            fprintf('CORE INVARIANTS: FAIL\n');
        end
        fprintf('------------------------------------------------------------\n\n');
    end

    if ~all(passed)
        failedNames = strjoin(testNames(~passed), '; ');
        error('test_core_invariants:Failed', ...
            '%d/%d core tests failed: %s', nnz(~passed), nTests, failedNames);
    end
end


function test_config_and_snr(o)
    cfg = build_fso_config(4, 1, 15, 0.1, ...
        'saMaxIter', o.saMaxIter, ...
        'saNStarts', o.saNStarts, ...
        'saUseParallel', false, ...
        'seedInit', o.seedInit, ...
        'logEvery', 0);

    expectedSigmaN2 = 1 / 10^(15/10);
    require_close(cfg.sigma_n_sq, expectedSigmaN2, 1e-14, ...
        'SNR convention mismatch.');
    require_true(cfg.M == 4, 'cfg.M mismatch.');
    require_close(cfg.P_avg, 1, 1e-14, 'cfg.P_avg mismatch.');
    require_close(cfg.sigma_X_sq, 0.1, 1e-14, 'cfg.sigma_X_sq mismatch.');
    require_true(isa(cfg.AMI_Evaluator,'function_handle'), 'Missing fast AMI evaluator.');
    require_true(isa(cfg.AMI_Validator,'function_handle'), 'Missing AMI validator.');
end


function test_uniform_pam(o)
    cfg = build_fso_config(8, 1, 20, 0.1, ...
        'saMaxIter', 1, 'saNStarts', o.saNStarts, ...
        'saUseParallel', false, 'seedInit', o.seedInit, 'logEvery', 0);

    x = define_constellation(8, 1, true, "mean");
    AMI_functions.assert_constellation_feasible(x, cfg, o.ProjectionTolerance);

    require_close(mean(x), 1, o.ProjectionTolerance, 'Uniform PAM mean power is not one.');
    require_true(all(x >= -o.ProjectionTolerance), 'Uniform PAM contains negative levels.');
    require_true(all(diff(x) >= cfg.SA.minGap-o.ProjectionTolerance), ...
        'Uniform PAM violates minGap.');
end


function test_projection(o)
    Ms = [4 8 16 32];
    s = RandStream('Threefry','Seed',424242);

    for M = Ms
        cfg = build_fso_config(M, 1, 15, 0.1, ...
            'saMaxIter', 1, 'saNStarts', o.saNStarts, ...
            'saUseParallel', false, 'minGap', 0.05, 'pinZero', true, ...
            'seedInit', o.seedInit, 'logEvery', 0);

        for trial = 1:o.ProjectionTrials
            raw = 3*randn(s,M,1) + 0.5*randn(s,1,1);
            x1 = AMI_functions.project_constellation_1D(raw, cfg);
            x2 = AMI_functions.project_constellation_1D(x1, cfg);

            AMI_functions.assert_constellation_feasible(x1, cfg, o.ProjectionTolerance);

            require_true(max(abs(x2-x1)) <= o.ProjectionTolerance, ...
                sprintf('Projection is not idempotent for M=%d, trial=%d.', M, trial));
            require_close(mean(x1), 1, o.ProjectionTolerance, ...
                sprintf('Mean-power constraint failed for M=%d.', M));
            require_true(abs(x1(1)) <= o.ProjectionTolerance, ...
                sprintf('pinZero failed for M=%d.', M));
            require_true(min(diff(x1)) >= cfg.SA.minGap-o.ProjectionTolerance, ...
                sprintf('minGap failed for M=%d.', M));
        end
    end
end


function test_projection_unpinned(o)
    cfg = build_fso_config(8, 1, 15, 0.1, ...
        'saMaxIter', 1, 'saNStarts', o.saNStarts, ...
        'saUseParallel', false, 'minGap', 0.02, 'pinZero', false, ...
        'seedInit', o.seedInit, 'logEvery', 0);

    raw = [-2.1; -0.4; 0.2; 0.21; 1.1; 2.8; 4.0; 5.5];
    x1 = AMI_functions.project_constellation_1D(raw, cfg);
    x2 = AMI_functions.project_constellation_1D(x1, cfg);

    AMI_functions.assert_constellation_feasible(x1, cfg, o.ProjectionTolerance);
    require_true(max(abs(x2-x1)) <= o.ProjectionTolerance, ...
        'Unpinned projection is not idempotent.');
    require_true(x1(1) >= -o.ProjectionTolerance, ...
        'Unpinned projection produced a negative first level.');
    require_close(mean(x1), 1, o.ProjectionTolerance, ...
        'Unpinned projection violates average power.');
end


function test_infeasible_gap()
    caught = false;
    try
        build_fso_config(32, 1, 15, 0.1, ...
            'minGap', 0.1, 'saMaxIter', 1, 'saNStarts', 2, ...
            'saUseParallel', false, 'seedInit', 1, 'logEvery', 0);
    catch ME
        caught = strcmp(ME.identifier, 'build_fso_config:InfeasibleMinGap');
        if ~caught
            rethrow(ME);
        end
    end
    require_true(caught, 'Infeasible minGap was not rejected.');
end


function test_ami_consistency(o)
    cfg = build_fso_config(4, 1, 15, 0.1, ...
        'saMaxIter', 1, 'saNStarts', o.saNStarts, ...
        'saUseParallel', false, 'seedInit', o.seedInit, 'logEvery', 0);
    x = define_constellation(4, 1, true, "mean");

    fast = cfg.AMI_Evaluator(x);
    val  = cfg.AMI_Validator(x);

    require_true(isfinite(fast) && isfinite(val), 'AMI evaluator returned a non-finite value.');
    require_true(fast >= -1e-12 && fast <= log2(cfg.M)+1e-8, 'Fast AMI outside theoretical range.');
    require_true(val  >= -1e-12 && val  <= log2(cfg.M)+1e-8, 'Validated AMI outside theoretical range.');
    require_true(abs(fast-val) <= o.AMITolerance, ...
        sprintf('|fast-validation|=%.3e exceeds tolerance %.3e.', abs(fast-val), o.AMITolerance));
end


function test_case_seed()
    s1 = fso_result_utils.case_seed(20260928, 4, 15, 0.1);
    s2 = fso_result_utils.case_seed(20260928, 4, 15, 0.1);
    s3 = fso_result_utils.case_seed(20260928, 4, 15, 0.2);
    s4 = fso_result_utils.case_seed(20260928, 8, 15, 0.1);

    require_true(s1 == s2, 'Identical channel cases produced different seeds.');
    require_true(s1 ~= s3, 'Different turbulence cases collided in case_seed().');
    require_true(s1 ~= s4, 'Different M cases collided in case_seed().');
end


function test_sa_repeatability(o)
    cfg = regression_sa_cfg(o);
    x0  = define_constellation(cfg.M, cfg.P_avg, true, "mean");

    [out1, r1] = sa_multistart(cfg, x0, '[REG-A]');
    [out2, r2] = sa_multistart(cfg, x0, '[REG-B]');

    tol = o.ReproTolerance;
    require_true(out1.masterSeed == out2.masterSeed, 'Same cfg produced different master seeds.');
    require_true(out1.bestStartFast == out2.bestStartFast, 'Fast winner is not reproducible.');
    require_true(out1.bestStartValidated == out2.bestStartValidated, ...
        'Validated winner is not reproducible.');
    require_close(out1.bestMIFast, out2.bestMIFast, tol, 'Fast best AMI is not reproducible.');
    require_close(out1.bestMIValidated, out2.bestMIValidated, tol, ...
        'Validated best AMI is not reproducible.');
    require_true(max(abs(out1.bestXValidated-out2.bestXValidated)) <= tol, ...
        'Validated winning constellation is not reproducible.');

    for k = 1:numel(r1)
        require_true(max(abs(r1(k).x0_raw-r2(k).x0_raw)) <= tol, ...
            sprintf('Restart %d initial point is not reproducible.', k));
        require_true(max(abs(r1(k).x_best-r2(k).x_best)) <= tol, ...
            sprintf('Restart %d best constellation is not reproducible.', k));
        require_close(r1(k).bestMIFast, r2(k).bestMIFast, tol, ...
            sprintf('Restart %d fast AMI is not reproducible.', k));
        require_close(r1(k).bestMIValidated, r2(k).bestMIValidated, tol, ...
            sprintf('Restart %d validated AMI is not reproducible.', k));

        AMI_functions.assert_constellation_feasible(r1(k).x_best, cfg, 1e-10);

        h = r1(k).history;
        require_true(~isempty(h.iter) && h.iter(1)==0 && h.iter(end)==cfg.SA.maxIter, ...
            sprintf('Restart %d history endpoints are invalid.', k));
        require_true(all(diff(h.iter) > 0), sprintf('Restart %d history iterations are not increasing.', k));
        require_true(all(diff(h.bestMI) >= -1e-12), ...
            sprintf('Restart %d best-MI history is not monotonic.', k));
        require_true(numel(h.iter)==numel(h.currentMI) && ...
                     numel(h.iter)==numel(h.bestMI) && ...
                     numel(h.iter)==numel(h.temperature) && ...
                     numel(h.iter)==numel(h.stepStd), ...
            sprintf('Restart %d history vectors have inconsistent lengths.', k));

        require_true(r1(k).nAccepted + r1(k).nRejected == cfg.SA.maxIter, ...
            sprintf('Restart %d acceptance accounting mismatch.', k));
        require_true(r1(k).nImprovingAccepted + r1(k).nWorseAccepted == r1(k).nAccepted, ...
            sprintf('Restart %d accepted-move accounting mismatch.', k));
        require_close(r1(k).acceptanceRate, r1(k).nAccepted/cfg.SA.maxIter, 1e-14, ...
            sprintf('Restart %d acceptance rate mismatch.', k));
    end
end


function test_validated_winner(o)
    cfg = regression_sa_cfg(o);
    x0  = define_constellation(cfg.M, cfg.P_avg, true, "mean");
    [out, results] = sa_multistart(cfg, x0, '[WINNER]');

    vals = [results.bestMIValidated];
    fast = [results.bestMIFast];
    [expectedVal, expectedIdx] = max(vals);
    [expectedFast, expectedFastIdx] = max(fast);

    require_true(out.bestStartValidated == expectedIdx, ...
        'Final restart index is not the maximum validated-AMI restart.');
    require_close(out.bestMIValidated, expectedVal, o.ReproTolerance, ...
        'Final validated AMI does not equal max restart validated AMI.');
    require_true(max(abs(out.bestXValidated-results(expectedIdx).x_best(:))) <= o.ReproTolerance, ...
        'Final validated constellation is not the validated-winning restart constellation.');

    require_true(out.bestStartFast == expectedFastIdx, 'Fast winner index mismatch.');
    require_close(out.bestMIFast, expectedFast, o.ReproTolerance, 'Fast winner AMI mismatch.');

    % Backward-compatible aliases must deliberately refer to validated winner.
    require_true(out.bestStart == out.bestStartValidated, 'out.bestStart alias is not validated winner.');
    require_close(out.bestMI, out.bestMIValidated, o.ReproTolerance, ...
        'out.bestMI alias is not validated AMI.');
    require_true(max(abs(out.bestX-out.bestXValidated)) <= o.ReproTolerance, ...
        'out.bestX alias is not validated constellation.');

    expectedChanged = expectedFastIdx ~= expectedIdx;
    require_true(out.selectionChanged == expectedChanged, 'selectionChanged flag is inconsistent.');
end


function cfg = regression_sa_cfg(o)
    cfg = build_fso_config(4, 1, 15, 0.1, ...
        'saMaxIter', o.saMaxIter, ...
        'saNStarts', o.saNStarts, ...
        'saUseParallel', false, ...
        'seedInit', o.seedInit, ...
        'logEvery', 0, ...
        'historyEvery', min(50,o.saMaxIter));
end


function require_true(condition, message)
    if ~condition
        error('test_core_invariants:AssertionFailed', '%s', message);
    end
end


function require_close(actual, expected, tol, message)
    if ~(isfinite(actual) && isfinite(expected) && abs(actual-expected) <= tol)
        error('test_core_invariants:AssertionFailed', ...
            '%s actual=%.16g expected=%.16g |delta|=%.3e tol=%.3e', ...
            message, actual, expected, abs(actual-expected), tol);
    end
end
