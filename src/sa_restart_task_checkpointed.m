function result = sa_restart_task_checkpointed(cfg,xBaseline,masterSeed,restartIndex,taskLabel,checkpointPath,signatureText,identity)
%SA_RESTART_TASK_CHECKPOINTED  Restart-task execution with atomic checkpointing.
%
% If checkpointPath exists and both the run signature and task identity match,
% the validated restart result is loaded and returned.  Otherwise the restart
% is executed through sa_restart_task(), then saved atomically via a temporary
% MAT file followed by movefile().  A partial/crashed task therefore never
% masquerades as a completed checkpoint.

    if nargin<6 || isempty(checkpointPath), checkpointPath=''; end
    if nargin<7, signatureText=''; end
    if nargin<8, identity=struct(); end
    checkpointPath=char(string(checkpointPath));
    signatureText=char(string(signatureText));

    if ~isempty(checkpointPath) && exist(checkpointPath,'file')==2
        try
            S=load(checkpointPath,'restartResult','checkpointMeta');
            ok=isfield(S,'restartResult') && isfield(S,'checkpointMeta') && ...
                isfield(S.checkpointMeta,'signatureText') && ...
                isfield(S.checkpointMeta,'identity') && ...
                strcmp(S.checkpointMeta.signatureText,signatureText) && ...
                isequaln(S.checkpointMeta.identity,identity);
            if ok
                result=S.restartResult;
                result.checkpointReused=true;
                result.checkpointPath=checkpointPath;
                return;
            end
        catch
            % Corrupt/incomplete/mismatched checkpoint: rerun and overwrite.
        end
    end

    result=sa_restart_task(cfg,xBaseline,masterSeed,restartIndex,taskLabel);

    if ~isempty(checkpointPath)
        checkpointMeta=struct();
        checkpointMeta.signatureText=signatureText;
        checkpointMeta.identity=identity;
        checkpointMeta.savedUTC=char(datetime('now','TimeZone','UTC', ...
            'Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
        restartResult=result; %#ok<NASGU>

        checkpointDir=fileparts(checkpointPath);
        if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end
        tmpPath=[tempname(checkpointDir) '.mat'];
        cleanup=onCleanup(@() cleanup_tmp(tmpPath)); %#ok<NASGU>
        save(tmpPath,'restartResult','checkpointMeta','-v7');
        [ok,msg]=movefile(tmpPath,checkpointPath,'f');
        if ~ok
            error('sa_restart_task_checkpointed:AtomicMoveFailed', ...
                'Could not finalize checkpoint %s: %s',checkpointPath,msg);
        end
    end

    result.checkpointReused=false;
    result.checkpointPath=checkpointPath;
end

function cleanup_tmp(path)
    if exist(path,'file')==2
        try, delete(path); catch, end
    end
end
