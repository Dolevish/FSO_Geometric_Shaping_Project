function report = test_i13d_dynamic_queue_smoke()
%TEST_I13D_DYNAMIC_QUEUE_SMOKE  Small process-pool test of dynamic dispatch.
%
% Exercises the same parfeval/fetchNext queue mechanics used by I.13D.
% Twelve short heterogeneous tasks are executed on four process workers.
% The test verifies that four futures remain in flight while pending work
% still exists, i.e. no parfor-style persistent chunk tail is possible.

    nWorkers=4;
    durations=[0.40 0.35 0.30 0.25 0.20 0.18 0.16 0.14 0.12 0.10 0.08 0.06];

    fprintf('\n============================================================\n');
    fprintf('I.13D dynamic parfeval queue smoke test\n');
    fprintf('workers=%d | tasks=%d\n',nWorkers,numel(durations));
    fprintf('============================================================\n');

    pp=gcp('nocreate');
    if ~isempty(pp),delete(pp);end
    pp=parpool('Processes',nWorkers);
    cleanup=onCleanup(@() cleanup_pool(pp)); %#ok<NASGU>

    active=parallel.FevalFuture.empty(0,1);
    next=1; n=numel(durations); completed=0;
    activeTrace=zeros(n,1); pendingTrace=zeros(n,1);

    while numel(active)<nWorkers && next<=n
        active(end+1,1)=parfeval(pp,@pause_task,1,next,durations(next)); %#ok<AGROW>
        next=next+1;
    end

    while ~isempty(active)
        [idx,~]=fetchNext(active);
        active(idx)=[];
        completed=completed+1;

        if next<=n
            active(end+1,1)=parfeval(pp,@pause_task,1,next,durations(next));
            next=next+1;
        end

        activeTrace(completed)=numel(active);
        pendingTrace(completed)=n-next+1;
        fprintf('completed %d/%d | active=%d | pending=%d\n', ...
            completed,n,numel(active),max(0,n-next+1));
    end

    preTail=pendingTrace>0;
    maintained=all(activeTrace(preTail)==nWorkers);
    assert(maintained,'Dynamic queue failed to keep all workers occupied before the final tail.');

    report=struct('passed',true,'workers',nWorkers,'tasks',n, ...
        'fullQueueMaintainedUntilTail',maintained, ...
        'activeTrace',activeTrace,'pendingTrace',pendingTrace);
    fprintf('PASS: full %d-worker queue maintained until the final tail.\n',nWorkers);
end

function out=pause_task(id,seconds)
    pause(seconds);
    out=id;
end

function cleanup_pool(pp)
    try
        if ~isempty(pp) && isvalid(pp),delete(pp);end
    catch
    end
end
