function report = test_k7_fig8_dual_export(varargin)
%TEST_K7_FIG8_DUAL_EXPORT Preflight: no optimization and no AMI evaluation.
%
% Verifies completed K6 data and exports BOTH K7 figures from synthetic
% data into an automatically removed temporary directory. Existing K6 MAT,
% checkpoints, tables and figures remain untouched.
%
% Usage:
%   test_k7_fig8_dual_export();
%   test_k7_fig8_dual_export('K6ResultsPath','/path/K6_repeatability_results.mat');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'K6ResultsPath','',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
    parse(p,varargin{:}); opts=p.Results;

    tasks={@()check_completed_k6(opts.K6ResultsPath), ...
        @check_synthetic_dual_figures, ...
        @check_invalid_statistics_rejected, ...
        @check_no_optimizer_dependency};
    labels={'Completed K6 source has 20 measured and consistent scores', ...
        'Separate figures: six files, zero-based main and zoom-only axes', ...
        'Sample-statistics mismatch is rejected', ...
        'Pure post-processing source; no SA/evaluator calls'};
    pass=false(numel(tasks),1);
    msg=strings(numel(tasks),1);

    fprintf('\n============================================================\n');
    fprintf('COMMIT K7: FIG. 8 DUAL-EXPORT PREFLIGHT (NO SA RUNS)\n');
    fprintf('============================================================\n');
    for j=1:numel(tasks)
        t=tic;
        try
            tasks{j}();
            pass(j)=true;msg(j)="PASS";
            fprintf('[%d/%d] %s ... PASS (%.2f s)\n', ...
                j,numel(tasks),labels{j},toc(t));
        catch ME
            msg(j)=string(ME.message);
            fprintf('[%d/%d] %s ... FAIL: %s\n', ...
                j,numel(tasks),labels{j},ME.message);
        end
    end
    report=struct('passed',pass,'messages',msg,'nPassed',nnz(pass), ...
        'nTotal',numel(pass));
    fprintf('\nK7 PREFLIGHT %d/%d PASS\n',nnz(pass),numel(pass));
    if ~all(pass)
        error('K7:PreflightFailed','K7 preflight failed; do not export yet.');
    end
end


function check_completed_k6(given)
    if strlength(string(given))>0
        path=char(string(given));
        if exist(path,'dir')==7
            path=fullfile(path,'K6_repeatability_results.mat');
        end
    else
        root=fullfile(fso_result_utils.revision2_results_root(), ...
            'commit_k6_fig8');
        d=dir(fullfile(root,'K6_*','K6_repeatability_results.mat'));
        if isempty(d)
            error('K7:MissingK6','No completed K6 MAT found.');
        end
        [~,i]=max([d.datenum]);
        path=fullfile(d(i).folder,d(i).name);
    end
    S=load(path,'R6');
    if ~isfield(S,'R6') || ...
            ~strcmp(string(S.R6.version),"K6-FIG8-I11p1-v1")
        error('K7:WrongK6','Wrong or missing K6 R6 struct.');
    end
    R=S.R6;
    if ~istable(R.restartTable)||height(R.restartTable)~=20 || ...
            ~isequal(R.restartTable.Restart(:),(1:20).') || ...
            ~iscell(R.restartResults)||numel(R.restartResults)~=20 || ...
            ~isstruct(R.statistics)||R.statistics.N~=20
        error('K7:IncompleteK6','K6 data not complete.');
    end
    y=R.restartTable.ValidatedAMI;
    for k=1:20
        if abs(R.restartResults{k}.bestMIValidated-y(k))>1e-10
            error('K7:K6Data','K6 MAT raw data/table mismatch.');
        end
    end
    if abs(R.statistics.MeanAMI-mean(y))>1e-10 || ...
            abs(R.statistics.SampleVarianceAMI-var(y,0))>1e-10 || ...
            abs(R.statistics.SampleStdAMI-std(y,0))>1e-10
        error('K7:SampleStats','K6 mean/variance/std mismatch.');
    end
    fprintf('  Source: %s\n',path);
end


function check_synthetic_dual_figures()
    [T,stats]=mock_data();
    Tbefore=T;Sbefore=stats;
    root=fullfile(fso_result_utils.revision2_results_root(), ...
        'commit_k7_fig8');
    if exist(root,'dir')~=7,mkdir(root);end
    tmp=tempname(root);mkdir(tmp);
    cleanup=onCleanup(@()cleanup_temp(tmp)); %#ok<NASGU>

    F=plot_k7_fig8_repeatability_dual(T,stats,tmp);
    figs=[F.main.figureHandle F.zoom.figureHandle];
    cleanupFigures=onCleanup(@()close_figs(figs)); %#ok<NASGU>

    required={'pdf','png','fig'};
    for tag={'main','zoom'}
        f=F.(tag{1});
        for j=1:numel(required)
            path=f.(required{j});
            if exist(path,'file')~=2
                error('K7:MissingArtifact','Missing %s',path);
            end
        end
        if numel(findall(f.figureHandle,'Type','axes'))~=1
            error('K7:ExtraAxes','Expected exactly one axes per figure, no inset.');
        end
        if f.axesHandle.Position(2)+f.axesHandle.Position(4)>0.79
            error('K7:StatsOcclude','Plot axes overlap statistics header.');
        end
    end
    if abs(F.main.yLimits(1))>eps || F.zoom.yLimits(1)<=0 || ...
            F.zoom.yLimits(1)>=min(T.ValidatedAMI)
        error('K7:AxisLimits','Unexpected main/zoom y-range.');
    end
    if ~isequaln(T,Tbefore)||~isequaln(stats,Sbefore)
        error('K7:InputMutation','Plot must not modify K6 input data.');
    end
end


function check_invalid_statistics_rejected()
    [T,S]=mock_data();
    S.SampleStdAMI=S.SampleStdAMI+0.05;
    failed=false;
    try
        plot_k7_fig8_repeatability_dual(T,S, ...
            fullfile(fso_result_utils.revision2_results_root(), ...
            'commit_k7_fig8','invalid_should_not_be_created'));
    catch ME
        failed=strcmp(ME.identifier,'K7:StatisticsMismatch');
    end
    if ~failed
        error('K7:MissingGuard','Incorrect sample statistics were not rejected.');
    end
end


function check_no_optimizer_dependency()
    src=fileread(which('export_k7_fig8_dual'));
    src2=fileread(which('plot_k7_fig8_repeatability_dual'));
    % Allow names of forbidden routines in COMMENTS, not callable tokens.
    lines=regexp(sprintf('%s\n%s',src,src2),'\r?\n','split');
    for k=1:numel(lines)
        line=strtrim(lines{k});
        if isempty(line)||startsWith(line,'%'),continue;end
        if ~isempty(regexp(line, ...
                '(^|[^A-Za-z_])(sa_restart_task|simulated_annealing|AMI_Evaluator|AMI_Validator)\s*\(', ...
                'once'))
            error('K7:OptimizerCall','Figure-only source calls optimization.');
        end
    end
end


function [T,S]=mock_data()
    y=1.704+0.00019*(1:20).' ...
        +0.0004*sin((1:20).'*1.83);
    T=table((1:20).',y,'VariableNames',{'Restart','ValidatedAMI'});
    S=struct('N',20,'MeanAMI',mean(y), ...
        'SampleVarianceAMI',var(y,0),'SampleStdAMI',std(y,0), ...
        'MinAMI',min(y),'MaxAMI',max(y), ...
        'RangeAMI',max(y)-min(y));
end
function close_figs(fs)
    for k=1:numel(fs),if isgraphics(fs(k)),close(fs(k));end,end
end
function cleanup_temp(folder)
    if exist(folder,'dir')==7,try,rmdir(folder,'s');catch,end,end
end
