function d = diagnose_i10p5a_m32_worker_scaling(varargin)
%DIAGNOSE_I10P5A_M32_WORKER_SCALING  Measure restart-task worker scaling.
%
% This diagnostic is intentionally short compared with production. It runs
% the same M=32 restart-task workload at several process-pool sizes so the
% final production worker count can be chosen from measured wall time rather
% than core count alone. M=32 is memory-bandwidth heavy, so more workers are
% not assumed to be faster.
%
% Default workload:
%   M=32, SNR=30 dB, sigma_X^2=0.3
%   8 restart tasks, 200 SA iterations each, one replicate
%   worker counts [2 4 6 8]
%
% The numerical results should be identical across worker counts because
% every restart is addressed by the same case seed and Threefry substream.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'WorkerCounts',[2 4 6 8],@(v) isnumeric(v)&&isvector(v)&&all(v>=1)&&all(mod(v,1)==0));
    addParameter(p,'MaxIter',200,@(v) isnumeric(v)&&isscalar(v)&&v>=1&&mod(v,1)==0);
    addParameter(p,'NStarts',8,@(v) isnumeric(v)&&isscalar(v)&&v>=1&&mod(v,1)==0);
    addParameter(p,'SigmaX2',0.3,@(v) isnumeric(v)&&isscalar(v)&&v>=0);
    addParameter(p,'RunLabel','I10p5a_worker_scaling',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:}); o=p.Results;

    wc=unique(double(o.WorkerCounts(:).'),'stable');
    rows=cell(numel(wc),1); referenceAMI=[]; referenceStarts=[];
    parent=fullfile(char(o.ResultsRoot),'i10p5a_worker_scaling');
    if ~exist(parent,'dir'),mkdir(parent);end

    for k=1:numel(wc)
        pool=gcp('nocreate');
        if ~isempty(pool)
            delete(pool);
        end
        parpool('Processes',wc(k));

        t=tic;
        b=sim_revision_AMI_vs_SNR_flat_restarts( ...
            'MVec',32,'SNRVec',30,'TurbulenceVec',double(o.SigmaX2), ...
            'saMaxIter',double(o.MaxIter),'saNStarts',double(o.NStarts),'nReplicates',1, ...
            'T0',0.4,'Tf',1e-3,'BaseStd0',0.15,'ItersPerTemp',50, ...
            'UseParallel',true,'ParallelWorkers',wc(k),'LongestFirst',true, ...
            'WriteCSV',false,'ResultsRoot',parent, ...
            'RunLabel',sprintf('%s_w%d',char(string(o.RunLabel)),wc(k)));
        wall=toc(t);

        s=load(b.caseFiles{1,1,1,1},'c'); starts=s.c.starts;
        ami=[starts.bestMIValidated];
        if isempty(referenceAMI)
            referenceAMI=s.c.winner.amiValidated;
            referenceStarts=ami;
        else
            assert(abs(s.c.winner.amiValidated-referenceAMI)<5e-12, ...
                'Validated winner changed with worker count.');
            assert(max(abs(ami-referenceStarts))<5e-12, ...
                'Per-restart validated AMI changed with worker count.');
        end

        rows{k}=struct('Workers',wc(k),'WallSeconds',wall, ...
            'SchedulerSeconds',b.parallel.schedulerWallSeconds, ...
            'WinnerAMI',s.c.winner.amiValidated, ...
            'MeanRestartSeconds',mean([starts.taskWallSeconds]), ...
            'MaxRestartSeconds',max([starts.taskWallSeconds]), ...
            'Tasks',b.parallel.nRestartTasks);
    end

    T=struct2table(vertcat(rows{:}));
    best=min(T.WallSeconds);
    T.SpeedupVsSlowest=max(T.WallSeconds)./T.WallSeconds;
    T.RelativeToBest=T.WallSeconds./best;

    d=struct(); d.results=T;
    [~,iBest]=min(T.WallSeconds);
    d.bestWorkers=T.Workers(iBest);
    d.bestWallSeconds=T.WallSeconds(iBest);
    d.note='Measured diagnostic only; do not call the selected worker count theoretically optimal.';

    out=fullfile(parent,[char(string(o.RunLabel)) '_summary.csv']);
    writetable(T,out);
    fprintf('\nI.10.5A M32 WORKER SCALING\n');
    disp(T);
    fprintf('Fastest tested worker count: %d | wall %.1fs\n',d.bestWorkers,d.bestWallSeconds);
end
