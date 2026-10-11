function out = export_k7_fig8_dual(source,varargin)
%EXPORT_K7_FIG8_DUAL Regenerate Fig. 8 from finished K6 WITHOUT rerunning SA.
%
% Input (any one):
%   out=export_k7_fig8_dual(R6);
%   out=export_k7_fig8_dual('/path/K6_repeatability_results.mat');
%   out=export_k7_fig8_dual();   % discovers latest completed K6 MAT
%
% Creates two completely separate 20-bar figures:
%   Fig8_SA_repeatability_main.{pdf,png,fig}  (zero-origin, no inset)
%   Fig8_SA_repeatability_zoom.{pdf,png,fig}  (magnified y-range)
%
% Does NOT call simulated_annealing, SA restart tasks, AMI_Evaluator, or
% AMI_Validator. Does NOT edit K6's original .mat, .csv or figures.
% New files go to results/ieee_revision2/commit_k7_fig8/K7_from_K6_<id>.

    if nargin<1,source=[];end
    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'CloseFigures',true,@logical_scalar);
    parse(p,varargin{:}); opt=p.Results;

    sourcePath='';
    if isempty(source)
        sourcePath=find_recent_k6_result();
        R=load_k6(sourcePath);
    elseif ischar(source)||(isstring(source)&&isscalar(source))
        sourcePath=char(string(source));
        if exist(sourcePath,'dir')==7
            sourcePath=fullfile(sourcePath,'K6_repeatability_results.mat');
        end
        R=load_k6(sourcePath);
    elseif isstruct(source)
        R=source;
        if isfield(R,'outputDirectory') && ...
                exist(fullfile(R.outputDirectory,'K6_repeatability_results.mat'),'file')==2
            sourcePath=fullfile(R.outputDirectory,'K6_repeatability_results.mat');
        end
    else
        error('K7:Input','Supply a completed K6 R6 struct or MAT file path.');
    end

    verify_k6(R,sourcePath);
    base=char(string(R.outputDirectory));
    [~,k6folder]=fileparts(base);
    if ~startsWith(k6folder,'K6_')
        error('K7:K6Directory','Expected K6_<UTC> source run folder.');
    end

    root=fso_result_utils.revision2_results_root();
    if strlength(string(opt.OutputDirectory))>0
        outputDir=char(string(opt.OutputDirectory));
    else
        outputDir=fullfile(root,'commit_k7_fig8', ...
            ['K7_from_' k6folder]);
    end
    target=canonical(outputDir);
    k7root=canonical(fullfile(root,'commit_k7_fig8'));
    if ~(strcmp(target,k7root)||startsWith(target,[k7root filesep]))
        error('K7:OutputDirectory', ...
            'K7 figures must be written under results/ieee_revision2/commit_k7_fig8.');
    end
    if exist(outputDir,'dir')~=7
        [ok,msg]=mkdir(outputDir);
        if ~ok,error('K7:mkdir','%s',msg);end
    end

    fprintf('\n============================================================\n');
    fprintf('COMMIT K7 - FIG. 8 DUAL EXPORT (READ-ONLY K6)\n');
    fprintf('Source K6: %s\n',base);
    fprintf('Confirmed 20 validated restart rows; 0 SA iterations scheduled.\n');
    fprintf('Mean %.9f | sample std %.9g | sample variance %.9g\n', ...
        R.statistics.MeanAMI, ...
        R.statistics.SampleStdAMI, ...
        R.statistics.SampleVarianceAMI);
    fprintf('Output K7: %s\n',outputDir);
    fprintf('============================================================\n');

    figureResult=plot_k7_fig8_repeatability_dual( ...
        R.restartTable,R.statistics,outputDir, ...
        'CloseFigures',opt.CloseFigures);

    out=struct();
    out.version='K7-FIG8-EXPORT-v1';
    out.sourceK6Directory=base;
    out.sourceK6Mat=sourcePath;
    out.sourceProduction=R.sourceProduction;
    out.outputDirectory=outputDir;
    out.nRestarts=20;
    out.scientificComputation='none; reuse 20 K6 validated AMI scores';
    out.statistics=R.statistics;
    out.main=figureResult.main;
    out.zoom=figureResult.zoom;

    % Do not persist figure handles (keep persistent artifact small).
    manifest=out;
    for tag={'main','zoom'}
        f=tag{1};
        h={'figureHandle','axesHandle'};
        present=intersect(h,fieldnames(manifest.(f)));
        if ~isempty(present)
            manifest.(f)=rmfield(manifest.(f),present);
        end
    end
    save(fullfile(outputDir,'K7_figure_manifest.mat'),'manifest','-v7');
    fprintf('K7 complete. Six figure files + manifest saved. K6 unchanged.\n');
end


function R=load_k6(path)
    if exist(path,'file')~=2
        error('K7:MissingK6','Completed K6 MAT not found: %s',path);
    end
    S=load(path,'R6');
    if ~isfield(S,'R6')
        error('K7:MissingR6','Expected variable R6 in: %s',path);
    end
    R=S.R6;
end


function path=find_recent_k6_result()
    root=fullfile(fso_result_utils.revision2_results_root(), ...
        'commit_k6_fig8');
    d=dir(fullfile(root,'K6_*','K6_repeatability_results.mat'));
    if isempty(d)
        error('K7:K6NotFound', ...
            'No finished K6 MAT found. Supply the path explicitly.');
    end
    [~,order]=sort([d.datenum],'descend');
    d=d(order);
    path=fullfile(d(1).folder,d(1).name);
    fprintf('K7 automatically selected latest completed K6: %s\n',path);
end


function verify_k6(R,path)
    required={'version','restartTable','statistics','restartResults', ...
        'outputDirectory','sourceProduction','signature'};
    if ~isstruct(R)||~all(isfield(R,required)) || ...
            ~strcmp(char(string(R.version)),'K6-FIG8-I11p1-v1')
        error('K7:WrongSource','Expected complete K6-FIG8-I11p1-v1 data.');
    end
    sig=R.signature;
    if ~isfield(sig,'M')||~isfield(sig,'SNRdB')|| ...
            ~isfield(sig,'SigmaR2')||~isfield(sig,'nRestarts')|| ...
            sig.M~=8 || sig.SNRdB~=20 || ...
            abs(sig.SigmaR2-.1)>1e-12 || sig.nRestarts~=20
        error('K7:ScientificPoint','K6 data are from the wrong operating point.');
    end
    T=R.restartTable;
    if ~istable(T)||height(T)~=20 || ...
            ~isequal(T.Restart(:),(1:20).') || ...
            ~ismember('ValidatedAMI',T.Properties.VariableNames)
        error('K7:K6RestartTable','Expected 20 ordered K6 restart scores.');
    end
    rr=R.restartResults;
    if ~iscell(rr)||numel(rr)~=20
        error('K7:K6RestartResults','Expected 20 original restart structs.');
    end
    for k=1:20
        v=rr{k};
        if ~isstruct(v)||~isfield(v,'bestMIValidated') || ...
                ~isfield(v,'startIndex') || ...
                v.startIndex~=k || ...
                ~isfinite(v.bestMIValidated) || ...
                abs(v.bestMIValidated-T.ValidatedAMI(k))>1e-10
            error('K7:K6DataMismatch', ...
                'Restart %d differs between K6 table and saved raw result.',k);
        end
    end
    st=R.statistics;
    if st.N~=20 || abs(mean(T.ValidatedAMI)-st.MeanAMI)>1e-10 || ...
            abs(var(T.ValidatedAMI,0)-st.SampleVarianceAMI)>1e-10 || ...
            abs(std(T.ValidatedAMI,0)-st.SampleStdAMI)>1e-10 || ...
            abs((max(T.ValidatedAMI)-min(T.ValidatedAMI))-st.RangeAMI)>1e-10
        error('K7:K6StatsMismatch','K6 raw scores and sample summary disagree.');
    end
    root=canonical(fullfile( ...
        fso_result_utils.revision2_results_root(),'commit_k6_fig8'));
    base=canonical(R.outputDirectory);
    if ~startsWith(base,[root filesep])
        error('K7:WrongSourceDirectory', ...
            'Input K6 run must reside under the original commit_k6_fig8 root.');
    end
    if ~isempty(path)&&~strcmp(canonical(fileparts(path)),base)
        error('K7:SourcePathMismatch', ...
            'Loaded MAT directory differs from R6.outputDirectory.');
    end
end


function x=canonical(path)
    x=char(java.io.File(char(string(path))).getCanonicalPath());
end
function tf=text_scalar(x)
    tf=ischar(x)||(isstring(x)&&isscalar(x));
end
function tf=logical_scalar(x)
    tf=islogical(x)&&isscalar(x);
end
