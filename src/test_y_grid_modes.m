function report = test_y_grid_modes(varargin)
%TEST_Y_GRID_MODES  Commit-I.6 regression suite for fixed/adaptive y-grid modes.
%
%   report = test_y_grid_modes()
%
% Tests:
%   1) fixed-bound remains the canonical default;
%   2) explicit fixed-bound reproduces the default evaluator;
%   3) candidate-adaptive preserves the exact constellation and does not
%      impose a hidden peak constraint;
%   4) representative adaptive AMI agrees with the independent validator;
%   5) repeated adaptive evaluation is deterministic;
%   6) metadata records yGridMode/yGridScale and GH rejects adaptive mode.
%
% Name-value options:
%   'AMITolerance'   max |adaptive-validator| (default 1e-3)
%   'Verbose'        print PASS/FAIL lines (default true)

    p=inputParser;
    addParameter(p,'AMITolerance',1e-3,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0);
    addParameter(p,'Verbose',true,@(v) islogical(v)&&isscalar(v));
    parse(p,varargin{:});
    o=p.Results;

    names={ ...
        'Canonical default remains fixed-bound'
        'Explicit fixed-bound reproduces default'
        'Adaptive mode does not clip or constrain xMax'
        'Adaptive fast-vs-validator accuracy'
        'Adaptive evaluator is deterministic'
        'Metadata and invalid combinations are explicit'};

    fns={ ...
        @test_default_mode
        @test_fixed_reproduction
        @test_no_hidden_peak_constraint
        @() test_accuracy(o.AMITolerance)
        @test_determinism
        @test_metadata_and_rejection};

    n=numel(names); passed=false(n,1); elapsed=nan(n,1); messages=strings(n,1);
    if o.Verbose
        fprintf('\n============================================================\n');
        fprintf('Commit I.6 - y-grid mode regression suite\n');
        fprintf('Adaptive accuracy tolerance = %.3e bit/symbol\n',o.AMITolerance);
        fprintf('============================================================\n');
    end

    for k=1:n
        if o.Verbose, fprintf('[%02d/%02d] %s ... ',k,n,names{k}); end
        t=tic;
        try
            fns{k}();
            elapsed(k)=toc(t); passed(k)=true; messages(k)="PASS";
            if o.Verbose, fprintf('PASS (%.2fs)\n',elapsed(k)); end
        catch ME
            elapsed(k)=toc(t); messages(k)=string(ME.message);
            if o.Verbose
                fprintf('FAIL (%.2fs)\n',elapsed(k));
                fprintf('        %s [%s]\n',ME.message,ME.identifier);
            end
        end
    end

    report=table(string(names),passed,elapsed,messages, ...
        'VariableNames',{'Test','Passed','ElapsedSeconds','Message'});

    if o.Verbose
        fprintf('\nPassed %d/%d tests in %.2fs.\n\n',nnz(passed),n,sum(elapsed,'omitnan'));
    end
    if ~all(passed)
        error('test_y_grid_modes:Failed','%d/%d y-grid tests failed.',nnz(~passed),n);
    end
end

function test_default_mode()
    cfg=build_fso_config(8,1,20,0.1,'minGap',0.01,'saMaxIter',1,'saNStarts',1);
    require_true(strcmp(cfg.yGridMode,'fixed-bound'),'Default yGridMode changed unexpectedly.');
    require_close(cfg.yGridScale,1.1,1e-14,'Default yGridScale mismatch.');
end

function test_fixed_reproduction()
    args={'minGap',0.01,'xMaxBound',7.79,'saMaxIter',1,'saNStarts',1, ...
        'FastFadingMethod','logh','loghDt',0.01,'loghSpanSigma',8};
    a=build_fso_config(8,1,20,0.1,args{:});
    b=build_fso_config(8,1,20,0.1,args{:},'YGridMode','fixed-bound','YGridScale',1.1);
    x=define_constellation(8,1,true,"mean");
    require_close(a.AMI_Evaluator(x),b.AMI_Evaluator(x),1e-13, ...
        'Explicit fixed-bound does not reproduce canonical default.');
end

function test_no_hidden_peak_constraint()
    % Feasible mean-one M=4 constellation with xMax well above xMaxBound.
    x=[0;0.01;0.02;3.97];
    cfg=build_fso_config(4,1,30,0.1, ...
        'minGap',0.01,'pinZero',true,'xMaxBound',2, ...
        'FastFadingMethod','logh','YGridMode','candidate-adaptive','YGridScale',1.1, ...
        'saMaxIter',1,'saNStarts',1);
    AMI_functions.assert_constellation_feasible(x,cfg,1e-10);
    before=x;
    v=cfg.AMI_Evaluator(x);
    require_true(isfinite(v),'Adaptive evaluator returned non-finite AMI.');
    require_true(isequal(x,before),'Adaptive evaluator modified the constellation.');
    require_true(max(x)>cfg.xMaxBound, ...
        'Test did not exercise xMax above the historical fixed bound.');
end

function test_accuracy(tol)
    cases=[8 20; 8 30; 32 20; 32 30];
    for k=1:size(cases,1)
        M=cases(k,1); snr=cases(k,2);
        xmax=fso_result_utils.feasible_xmax_bound(M,1,0.01,true);
        cfg=build_fso_config(M,1,snr,0.1, ...
            'minGap',0.01,'pinZero',true,'xMaxBound',xmax, ...
            'FastFadingMethod','logh','loghDt',0.01,'loghSpanSigma',8, ...
            'YGridMode','candidate-adaptive','YGridScale',1.1, ...
            'saMaxIter',1,'saNStarts',1);
        x=define_constellation(M,1,true,"mean");
        fast=cfg.AMI_Evaluator(x);
        val=cfg.AMI_Validator(x);
        require_true(abs(fast-val)<=tol, ...
            sprintf('M=%d SNR=%g: |fast-validator|=%.3e > %.3e.',M,snr,abs(fast-val),tol));
    end
end

function test_determinism()
    M=16; xmax=fso_result_utils.feasible_xmax_bound(M,1,0.01,true);
    cfg=build_fso_config(M,1,30,0.1, ...
        'minGap',0.01,'xMaxBound',xmax,'FastFadingMethod','logh', ...
        'YGridMode','candidate-adaptive','YGridScale',1.1, ...
        'saMaxIter',1,'saNStarts',1);
    raw=linspace(0,4,M).';
    x=AMI_functions.project_constellation_1D(raw,cfg);
    a=cfg.AMI_Evaluator(x); b=cfg.AMI_Evaluator(x);
    require_close(a,b,1e-13,'Repeated adaptive evaluation is not deterministic.');
end

function test_metadata_and_rejection()
    cfg=build_fso_config(8,1,20,0.1, ...
        'minGap',0.01,'FastFadingMethod','logh', ...
        'YGridMode','candidate-adaptive','YGridScale',1.1, ...
        'saMaxIter',1,'saNStarts',1);
    snap=fso_result_utils.config_snapshot(cfg);
    require_true(isfield(snap,'yGridMode')&&strcmp(snap.yGridMode,'candidate-adaptive'), ...
        'Config snapshot does not record adaptive yGridMode.');
    require_true(isfield(snap,'yGridScale'), 'Config snapshot does not record yGridScale.');
    require_close(snap.yGridScale,1.1,1e-14,'Snapshot yGridScale mismatch.');

    caught=false;
    try
        build_fso_config(8,1,20,0.1,'FastFadingMethod','gh', ...
            'YGridMode','candidate-adaptive','saMaxIter',1,'saNStarts',1);
    catch ME
        caught=strcmp(ME.identifier,'build_fso_config:AdaptiveYGridRequiresLogh');
        if ~caught, rethrow(ME); end
    end
    require_true(caught,'GH + candidate-adaptive mode was not rejected explicitly.');
end

function require_true(tf,msg)
    if ~tf, error('test_y_grid_modes:RequirementFailed','%s',msg); end
end

function require_close(a,b,tol,msg)
    if ~(isfinite(a)&&isfinite(b)&&abs(a-b)<=tol)
        error('test_y_grid_modes:RequirementFailed','%s | a=%.16g b=%.16g diff=%.3e', ...
            msg,a,b,abs(a-b));
    end
end
