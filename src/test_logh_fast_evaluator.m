function report = test_logh_fast_evaluator(varargin)
%TEST_LOGH_FAST_EVALUATOR  Commit-I.3 regression tests for fast log-h AMI.
%
%   report = test_logh_fast_evaluator()
%
% The test is intentionally independent of SA.  It verifies:
%   1) build_fso_config defaults to the production log-h evaluator;
%   2) log-h fast AMI agrees with the independent validator on representative
%      low/mid/high-SNR Uniform-PAM cases;
%   3) the no-turbulence branch is consistent;
%   4) explicit legacy GH selection remains available and matches a direct
%      call to AMI_functions.AMI_noCSI_fast_grid();
%   5) config snapshots preserve the fast-integration metadata;
%   6) invalid fast-method names are rejected.
%
% Name-value options:
%   'AMITolerance'     allowed |fast-validator| (default 1e-3)
%   'DefaultDt'        expected production log-h dt (default 0.01)
%   'DefaultSpanSigma' expected production span (default 8)
%   'Verbose'          print PASS/FAIL details (default true)
%
% This test does not write result files and does not modify RNG state.

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p,'AMITolerance',1e-3,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p,'DefaultDt',0.01,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p,'DefaultSpanSigma',8,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p,'Verbose',true,@(v) islogical(v) && isscalar(v));
    parse(p,varargin{:});
    o = p.Results;

    testNames = {
        'Canonical default is uniform log-h'
        'Representative fast-vs-validator accuracy'
        'No-turbulence consistency'
        'Legacy GH path remains selectable'
        'Config snapshot records integration metadata'
        'Invalid fast method is rejected'
        };

    testFns = {
        @() test_default(o)
        @() test_accuracy(o)
        @() test_awgn(o)
        @() test_legacy_gh()
        @() test_snapshot(o)
        @() test_bad_method()
        };

    n = numel(testNames);
    passed = false(n,1);
    elapsed = nan(n,1);
    messages = strings(n,1);

    if o.Verbose
        fprintf('\n============================================================\n');
        fprintf('Commit I.3 - fast log-h AMI regression suite\n');
        fprintf('AMI tolerance=%.3e | dt=%.5g | span=+/-%.3g sigma_t\n', ...
            o.AMITolerance,o.DefaultDt,o.DefaultSpanSigma);
        fprintf('============================================================\n');
    end

    for k=1:n
        if o.Verbose, fprintf('[%02d/%02d] %s ... ',k,n,testNames{k}); end
        t=tic;
        try
            testFns{k}();
            elapsed(k)=toc(t);
            passed(k)=true;
            messages(k)="PASS";
            if o.Verbose, fprintf('PASS (%.2fs)\n',elapsed(k)); end
        catch ME
            elapsed(k)=toc(t);
            messages(k)=string(ME.message);
            if o.Verbose
                fprintf('FAIL (%.2fs)\n',elapsed(k));
                fprintf('        %s [%s]\n',ME.message,ME.identifier);
            end
        end
    end

    report=table(string(testNames),passed,elapsed,messages, ...
        'VariableNames',{'Test','Passed','ElapsedSeconds','Message'});

    if o.Verbose
        fprintf('\nPassed %d/%d tests in %.2fs.\n',nnz(passed),n,sum(elapsed,'omitnan'));
    end
    if ~all(passed)
        error('test_logh_fast_evaluator:Failed','%d/%d Commit-I.3 tests failed.',nnz(~passed),n);
    end
end


function test_default(o)
    cfg=build_fso_config(8,1,20,0.1, ...
        'minGap',0.01,'saMaxIter',1,'saNStarts',1,'logEvery',0);
    require_true(strcmp(cfg.fastFadingMethod,'logh'),'Canonical fast method is not logh.');
    require_close(cfg.loghDt,o.DefaultDt,1e-14,'Default loghDt mismatch.');
    require_close(cfg.loghSpanSigma,o.DefaultSpanSigma,1e-14,'Default loghSpanSigma mismatch.');
    require_true(isa(cfg.AMI_Evaluator,'function_handle'),'Missing fast evaluator handle.');
end


function test_accuracy(o)
    % Covers low, mid and high SNR including the difficult M=32 high-SNR case.
    cases=[ ...
         8   5   0.1
         8  20   0.1
         8  30   0.1
        32  20   0.1
        32  30   0.1];

    for k=1:size(cases,1)
        M=cases(k,1); snr=cases(k,2); sig=cases(k,3);
        xmax=fso_result_utils.feasible_xmax_bound(M,1,0.01,true);
        cfg=build_fso_config(M,1,snr,sig, ...
            'FastFadingMethod','logh','loghDt',o.DefaultDt, ...
            'loghSpanSigma',o.DefaultSpanSigma,'xMaxBound',xmax, ...
            'minGap',0.01,'saMaxIter',1,'saNStarts',1,'logEvery',0);
        x=define_constellation(M,1,true,"mean");
        fast=cfg.AMI_Evaluator(x);
        val=cfg.AMI_Validator(x);
        require_true(isfinite(fast) && isfinite(val),'Non-finite AMI value.');
        require_true(abs(fast-val)<=o.AMITolerance, ...
            sprintf('M=%d SNR=%.1f: |fast-val|=%.3e > %.3e.', ...
            M,snr,abs(fast-val),o.AMITolerance));
    end
end


function test_awgn(o)
    M=8; xmax=fso_result_utils.feasible_xmax_bound(M,1,0.01,true);
    common={'xMaxBound',xmax,'minGap',0.01,'saMaxIter',1,'saNStarts',1,'logEvery',0};
    cfgL=build_fso_config(M,1,20,0,common{:}, ...
        'FastFadingMethod','logh','loghDt',o.DefaultDt,'loghSpanSigma',o.DefaultSpanSigma);
    cfgG=build_fso_config(M,1,20,0,common{:}, ...
        'FastFadingMethod','gh','ghN_h',40);
    x=define_constellation(M,1,true,"mean");
    a=cfgL.AMI_Evaluator(x);
    b=cfgG.AMI_Evaluator(x);
    require_close(a,b,1e-12,'AWGN log-h and GH fast paths disagree.');
    val=cfgL.AMI_Validator(x);
    require_true(abs(a-val)<=o.AMITolerance,'AWGN fast-validator mismatch exceeds tolerance.');
end


function test_legacy_gh()
    M=8; xmax=fso_result_utils.feasible_xmax_bound(M,1,0.01,true);
    cfg=build_fso_config(M,1,20,0.1, ...
        'FastFadingMethod','gh','ghN_h',80,'xMaxBound',xmax, ...
        'minGap',0.01,'saMaxIter',1,'saNStarts',1,'logEvery',0);
    require_true(strcmp(cfg.fastFadingMethod,'gh'),'Explicit GH selection was not preserved.');
    x=define_constellation(M,1,true,"mean");
    viaCfg=cfg.AMI_Evaluator(x);
    direct=AMI_functions.AMI_noCSI_fast_grid(x,cfg.px,cfg,80,cfg.y_grid);
    require_close(viaCfg,direct,1e-12,'Configured GH handle differs from direct GH evaluator.');
end


function test_snapshot(o)
    cfg=build_fso_config(8,1,20,0.1, ...
        'FastFadingMethod','logh','loghDt',o.DefaultDt, ...
        'loghSpanSigma',o.DefaultSpanSigma,'loghYBlockSize',512, ...
        'minGap',0.01,'saMaxIter',1,'saNStarts',1,'logEvery',0);
    snap=fso_result_utils.config_snapshot(cfg);
    require_true(strcmp(snap.fastFadingMethod,'logh'),'Snapshot lost fastFadingMethod.');
    require_close(snap.loghDt,o.DefaultDt,1e-14,'Snapshot lost loghDt.');
    require_close(snap.loghSpanSigma,o.DefaultSpanSigma,1e-14,'Snapshot lost loghSpanSigma.');
    require_true(snap.loghYBlockSize==512,'Snapshot lost loghYBlockSize.');
    require_true(isfield(snap,'yGridInfo'),'Snapshot missing y-grid metadata.');
    require_true(~isfield(snap,'AMI_Evaluator') && ~isfield(snap,'AMI_Validator'), ...
        'Snapshot retained function handles.');
end


function test_bad_method()
    caught=false;
    try
        build_fso_config(4,1,15,0.1,'FastFadingMethod','not-a-method', ...
            'saMaxIter',1,'saNStarts',1,'logEvery',0);
    catch ME
        caught=strcmp(ME.identifier,'build_fso_config:BadFastFadingMethod');
        if ~caught, rethrow(ME); end
    end
    require_true(caught,'Invalid fast method was not rejected.');
end


function require_true(tf,msg)
    if ~tf, error('test_logh_fast_evaluator:AssertionFailed','%s',msg); end
end


function require_close(a,b,tol,msg)
    if ~(isfinite(a) && isfinite(b) && abs(a-b)<=tol)
        error('test_logh_fast_evaluator:AssertionFailed', ...
            '%s | a=%.16g b=%.16g |diff|=%.3e tol=%.3e',msg,a,b,abs(a-b),tol);
    end
end
