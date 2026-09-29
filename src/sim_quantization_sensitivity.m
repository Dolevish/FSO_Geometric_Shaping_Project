function qbundle = sim_quantization_sensitivity(source,varargin)
%SIM_QUANTIZATION_SENSITIVITY  Finite-resolution implementation study.
%
% Commit H: reviewer-response experiment for practical implementation and
% quantization sensitivity of optimized FSO geometric-shaping constellations.
%
%   qbundle = sim_quantization_sensitivity(source)
%   qbundle = sim_quantization_sensitivity(source,'Name',Value,...)
%
% SOURCE may be either:
%   1) a minimum-spacing bundle returned by sim_min_gap_sensitivity(), or
%   2) the path to a .mat file containing that bundle as variable 'bundle'.
%
% The experiment applies a uniform B-bit DAC quantizer to each independently
% optimized source constellation.  The DAC full-scale is
%
%       X_FS = FullScaleMargin * max(x)
%
% and each level is mapped to the nearest code in {0,...,2^B-1}.  By default
% the quantized constellation is then multiplied by one scalar gain so that
%
%       mean(x_q) = P_avg.
%
% This models finite digital amplitude resolution followed by calibrated
% optical gain.  The scalar calibration preserves DAC-code collisions and
% relative code geometry while ensuring a fair AMI comparison at the same
% average optical power.  Set PowerCalibrate=false only for diagnostics.
%
% Primary outputs per M / d_min / B / independent replicate:
%   * independently validated quantized GS AMI;
%   * signed AMI penalty and absolute AMI change relative to the SAME source
%     optimized constellation (paired quantization comparison);
%   * quantized shaping gain against equally quantized Uniform PAM;
%   * number of unique DAC levels and number of collisions;
%   * post-calibration minimum adjacent spacing and quantization step;
%   * RMS/max constellation displacement, peak level and gain calibration;
%   * pre/post-calibration average-power error and clipping count.
%
% Default bit depths:
%       B = [4 5 6 7 8 10]
%
% Name-value options:
%   'MVec'             source M values to analyze ([] -> all in source)
%   'MinGapVec'        source d_min values ([] -> all in source)
%   'ReplicateVec'     source replicate indices ([] -> all in source)
%   'BitsVec'          DAC resolutions (default [4 5 6 7 8 10])
%   'FullScaleMargin'  X_FS/max(x), >=1 recommended (default 1)
%   'PowerCalibrate'   restore mean power after quantization (default true)
%   'LossThresholds'   abs-AMI thresholds used in summary (default [1e-2 1e-3])
%   'UseParallelCases' parallelize source-constellation cases (default true)
%   'MakePlots'        create summary figures (default true)
%   'SaveFigures'      save PNG/FIG when plotting (default true)
%   'RunLabel'         optional run label (default '')
%   'ResultsRoot'      default repository/results/ieee_revision
%
% Methodological notes:
%   * Quantization is paired within each source replicate.  Therefore AMI
%     change is not contaminated by optimizer-to-optimizer variability.
%   * All reported GS AMIs use cfg.AMI_Validator (independent validator), not
%     the fast GH objective.
%   * Uniform PAM is quantized using the same bit depth, full-scale rule and
%     power-calibration rule before computing quantized shaping gain.
%   * d_min=0 source constellations may contain nearly coincident levels; the
%     unique-level metric explicitly reveals when finite DAC resolution turns
%     those near-collisions into exact code collisions.
%
% This script does not modify or rerun the SA optimizer.  It consumes saved
% optimized constellations from the Commit-F minimum-spacing study.

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p,'source',@(v) isstruct(v) || ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'MVec',[],@(v) isempty(v) || (isnumeric(v) && isvector(v) && all(isfinite(v))));
    addParameter(p,'MinGapVec',[],@(v) isempty(v) || (isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>=0)));
    addParameter(p,'ReplicateVec',[],@(v) isempty(v) || (isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>=1) && all(mod(v,1)==0)));
    addParameter(p,'BitsVec',[4 5 6 7 8 10],@(v) isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>=1) && all(v<=16) && all(mod(v,1)==0));
    addParameter(p,'FullScaleMargin',1,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
    addParameter(p,'PowerCalibrate',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'LossThresholds',[1e-2 1e-3],@(v) isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v>0));
    addParameter(p,'UseParallelCases',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'MakePlots',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'SaveFigures',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v) || isstring(v));
    parse(p,source,varargin{:});
    o = p.Results;

    [src,sourceDescription] = load_source_bundle(source);
    validate_source_bundle(src);

    srcM   = double(src.axes.M(:).');
    srcGap = double(src.axes.minGap(:).');
    srcRep = double(src.axes.replicate(:).');

    if isempty(o.MVec), MVec = srcM; else, MVec = unique(double(o.MVec(:).'),'stable'); end
    if isempty(o.MinGapVec), gapVec = srcGap; else, gapVec = unique(double(o.MinGapVec(:).'),'stable'); end
    if isempty(o.ReplicateVec), repVec = srcRep; else, repVec = unique(double(o.ReplicateVec(:).'),'stable'); end
    bitsVec = unique(double(o.BitsVec(:).'),'stable');
    thresholds = unique(double(o.LossThresholds(:).'),'stable');

    idxM   = map_axis_values(MVec,srcM,'M');
    idxGap = map_axis_values(gapVec,srcGap,'d_min');
    idxRep = map_axis_values(repVec,srcRep,'replicate');

    nM = numel(MVec); nG = numel(gapVec); nB = numel(bitsVec); nR = numel(repVec);

    SNR_dB      = double(src.axes.SNR_dB);
    sigma_X_sq  = double(src.axes.sigma_X_sq);
    P_avg       = double(src.axes.P_avg);

    % ------------------------------------------------------------------
    % Immutable run identity / no-overwrite output directory.
    % ------------------------------------------------------------------
    resultsRoot = char(o.ResultsRoot);
    conditionTag = sprintf('SNR%s_sig%s',num_token(SNR_dB,2),num_token(sigma_X_sq,4));
    conditionDir = fullfile(resultsRoot,'quantization_sensitivity',conditionTag);

    runMeta = fso_result_utils.run_metadata(mfilename);
    labelToken = sanitize_run_label(o.RunLabel);
    bitsToken = sprintf('B%s',join_numeric_token(bitsVec));
    calToken = sprintf('cal%d',logical(o.PowerCalibrate));
    fsToken = sprintf('FS%s',num_token(o.FullScaleMargin,3));
    utcToken = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requestedRunId = sprintf('%s_%s_%s_%s_%s',labelToken,bitsToken,fsToken,calToken,utcToken);
    [outDir,runId] = make_unique_run_directory(conditionDir,requestedRunId);
    caseDir = fullfile(outDir,'cases');
    figDir = fullfile(outDir,'figures');
    mkdir(caseDir);
    if o.MakePlots && o.SaveFigures, mkdir(figDir); end

    runMeta.runId = runId;
    runMeta.runLabel = char(string(o.RunLabel));
    runMeta.outputDirectory = outDir;
    runMeta.conditionTag = conditionTag;

    % ------------------------------------------------------------------
    % Quantized Uniform-PAM baselines.  These depend only on M and B.
    % ------------------------------------------------------------------
    baseline = struct();
    baseline.xPAM = cell(nM,1);
    baseline.amiUnquantized = nan(nM,1);
    baseline.amiQuantized = nan(nM,nB);
    baseline.uniqueLevels = nan(nM,nB);
    baseline.collisionCount = nan(nM,nB);
    baseline.quantStep = nan(nM,nB);

    fprintf('\n============================================================\n');
    fprintf('IEEE revision quantization / implementation sensitivity\n');
    fprintf('Run ID: %s\n',runId);
    if strlength(string(o.RunLabel)) > 0
        fprintf('Run label: %s\n',char(string(o.RunLabel)));
    end
    fprintf('Source: %s\n',sourceDescription);
    fprintf('M = %s | d_min = %s | replicates = %s\n',mat2str(MVec),mat2str(gapVec),mat2str(repVec));
    fprintf('SNR=%.2f dB | sigma_X^2=%.4f | Pavg=%.4g\n',SNR_dB,sigma_X_sq,P_avg);
    fprintf('DAC bits = %s | full-scale margin=%.4g | power calibration=%d\n', ...
        mat2str(bitsVec),o.FullScaleMargin,o.PowerCalibrate);
    fprintf('Total source constellations: %d | quantized GS evaluations: %d\n', ...
        nM*nG*nR,nM*nG*nR*nB);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    for iM = 1:nM
        M = MVec(iM);
        cfg = validator_config(M,P_avg,SNR_dB,sigma_X_sq);
        xPAM = define_constellation(M,P_avg,true,"mean");
        baseline.xPAM{iM} = xPAM(:);
        baseline.amiUnquantized(iM) = cfg.AMI_Validator(xPAM(:));
        for iB = 1:nB
            [xq,qi] = quantize_dac(xPAM(:),bitsVec(iB),P_avg,o.FullScaleMargin,o.PowerCalibrate);
            baseline.amiQuantized(iM,iB) = cfg.AMI_Validator(xq);
            baseline.uniqueLevels(iM,iB) = qi.nUniqueLevels;
            baseline.collisionCount(iM,iB) = qi.collisionCount;
            baseline.quantStep(iM,iB) = qi.quantStepAfterCalibration;
        end
    end

    [iMGrid,iGGrid,iRGrid] = ndgrid(1:nM,1:nG,1:nR);
    taskM = iMGrid(:); taskG = iGGrid(:); taskR = iRGrid(:);
    nCases = numel(taskM);
    caseResults = cell(nCases,1);
    caseFiles = cell(nCases,1);

    doPar = logical(o.UseParallelCases);
    if doPar
        try
            doPar = license('test','Distrib_Computing_Toolbox') ~= 0;
        catch
            doPar = false;
        end
    end

    if doPar
        pool = gcp('nocreate');
        if isempty(pool), pool = parpool('Processes'); end
        fprintf('Outer case parallelism: %d process workers\n',pool.NumWorkers);
        Q = parallel.pool.DataQueue;
        doneCount = 0;
        afterEach(Q,@(~) progress_tick());

        parfor k = 1:nCases
            [caseResults{k},caseFiles{k}] = run_quant_case( ...
                taskM(k),taskG(k),taskR(k),MVec,gapVec,repVec,idxM,idxGap,idxRep, ...
                bitsVec,src,baseline,o,runMeta,caseDir,P_avg,SNR_dB,sigma_X_sq);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k = 1:nCases
            [caseResults{k},caseFiles{k}] = run_quant_case( ...
                taskM(k),taskG(k),taskR(k),MVec,gapVec,repVec,idxM,idxGap,idxRep, ...
                bitsVec,src,baseline,o,runMeta,caseDir,P_avg,SNR_dB,sigma_X_sq);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    caseGrid = reshape(caseResults,nM,nG,nR);
    fileGrid = reshape(caseFiles,nM,nG,nR);
    summary = build_summary(caseGrid,MVec,gapVec,bitsVec,repVec,thresholds);

    qbundle = struct();
    qbundle.meta = runMeta;
    qbundle.experiment = 'quantization_sensitivity';
    qbundle.run = struct('id',runId,'label',char(string(o.RunLabel)), ...
        'conditionTag',conditionTag,'outputDirectory',outDir);
    qbundle.source = struct('description',sourceDescription, ...
        'runId',source_run_id(src),'experiment',source_experiment(src));
    qbundle.axes = struct('M',MVec,'minGap',gapVec,'bits',bitsVec,'replicate',repVec, ...
        'SNR_dB',SNR_dB,'sigma_X_sq',sigma_X_sq,'P_avg',P_avg);
    qbundle.settings = struct('fullScaleMargin',double(o.FullScaleMargin), ...
        'powerCalibrate',logical(o.PowerCalibrate), ...
        'lossThresholds',thresholds,'pairedWithinSourceReplicate',true);
    qbundle.baseline = baseline;
    qbundle.caseFiles = fileGrid;
    qbundle.summary = summary;
    qbundle.figures = {};

    if o.MakePlots
        qbundle.figures = make_summary_plots(qbundle,figDir,logical(o.SaveFigures));
    end

    bundlePath = fullfile(outDir,'quantizationSensitivity.mat');
    save(bundlePath,'qbundle','-v7.3');
    fprintf('\nSaved quantization bundle: %s\n',bundlePath);
    print_summary(qbundle);

    function progress_tick()
        doneCount = doneCount + 1;
        fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function [c,path] = run_quant_case(iM,iG,iR,MVec,gapVec,repVec,idxM,idxGap,idxRep, ...
        bitsVec,src,baseline,o,runMeta,caseDir,P_avg,SNR_dB,sigma_X_sq)

    M = MVec(iM); d = gapVec(iG); repLabel = repVec(iR);
    sM = idxM(iM); sG = idxGap(iG); sR = idxRep(iR);

    x = src.summary.xWinner{sM,sG,sR};
    x = real(x(:));
    storedAMI = src.summary.amiGSValidated(sM,sG,sR);

    cfg = validator_config(M,P_avg,SNR_dB,sigma_X_sq);
    sourceAMI = cfg.AMI_Validator(x);

    nB = numel(bitsVec);
    q = repmat(struct(),nB,1);
    for iB = 1:nB
        B = bitsVec(iB);
        [xq,qi] = quantize_dac(x,B,P_avg,o.FullScaleMargin,o.PowerCalibrate);
        amiQ = cfg.AMI_Validator(xq);

        qi.bits = B;
        qi.x = xq;
        qi.amiValidated = amiQ;
        qi.amiPenalty = sourceAMI - amiQ;
        qi.deltaAMI = amiQ - sourceAMI;
        qi.absAMIChange = abs(amiQ-sourceAMI);
        qi.quantizedPAMAMI = baseline.amiQuantized(iM,iB);
        qi.quantizedShapingGain = amiQ - baseline.amiQuantized(iM,iB);
        q(iB) = qi;
    end

    c = struct();
    c.meta = runMeta;
    c.case = struct('M',M,'minGap',d,'replicate',repLabel, ...
        'SNR_dB',SNR_dB,'sigma_X_sq',sigma_X_sq,'P_avg',P_avg);
    c.source = struct('x',x,'storedAMIValidated',storedAMI, ...
        'revalidatedAMI',sourceAMI,'revalidationGap',sourceAMI-storedAMI, ...
        'actualMinGap',safe_min_diff(x),'xMax',max(x));
    c.quantized = q;
    c.baseline = struct('xPAM',baseline.xPAM{iM}, ...
        'amiPAMUnquantized',baseline.amiUnquantized(iM), ...
        'amiPAMQuantized',baseline.amiQuantized(iM,:));

    path = fullfile(caseDir,sprintf('case_M%d_gap%s_rep%02d.mat', ...
        M,num_token(d,4),repLabel));
    save(path,'c','-v7.3');
end


function [xq,info] = quantize_dac(x,B,P_avg,margin,powerCalibrate)
    x = real(x(:));
    if isempty(x) || any(~isfinite(x)) || any(x < -1e-12)
        error('sim_quantization_sensitivity:BadConstellation', ...
            'Source constellation must be finite and nonnegative.');
    end
    x = max(x,0);
    if max(x) <= 0
        error('sim_quantization_sensitivity:ZeroConstellation', ...
            'Cannot quantize an all-zero constellation.');
    end

    nCodes = 2^double(B);
    maxCode = nCodes-1;
    fullScale = double(margin)*max(x);

    u = x/fullScale;
    clippedLow = u < 0;
    clippedHigh = u > 1;
    u = min(max(u,0),1);
    codes = round(u*maxCode);
    xRaw = (codes/maxCode)*fullScale;

    meanBefore = mean(xRaw);
    if meanBefore <= eps
        error('sim_quantization_sensitivity:DegenerateQuantization', ...
            'Quantization collapsed all levels to zero.');
    end

    if powerCalibrate
        scale = P_avg/meanBefore;
    else
        scale = 1;
    end
    xq = xRaw*scale;

    uniqueCodes = unique(codes,'stable');
    nUnique = numel(uniqueCodes);
    uniqueX = unique(xq,'sorted');

    info = struct();
    info.codes = codes;
    info.fullScaleBeforeCalibration = fullScale;
    info.powerCalibrationScale = scale;
    info.quantStepBeforeCalibration = fullScale/maxCode;
    info.quantStepAfterCalibration = fullScale*scale/maxCode;
    info.meanPowerBeforeCalibration = meanBefore;
    info.meanPowerAfterCalibration = mean(xq);
    info.meanPowerErrorBeforeCalibration = meanBefore-P_avg;
    info.meanPowerErrorAfterCalibration = mean(xq)-P_avg;
    info.nUniqueLevels = nUnique;
    info.collisionCount = numel(x)-nUnique;
    info.allLevelsDistinct = (nUnique == numel(x));
    info.minAdjacentSpacing = safe_min_diff(xq);
    if numel(uniqueX) > 1
        info.minDistinctSpacing = min(diff(uniqueX));
    else
        info.minDistinctSpacing = 0;
    end
    info.xMax = max(xq);
    info.rmsDisplacement = sqrt(mean((xq-x).^2));
    info.maxAbsDisplacement = max(abs(xq-x));
    info.clippingCount = sum(clippedLow | clippedHigh);
end


function cfg = validator_config(M,P_avg,SNR_dB,sigma_X_sq)
    % GH settings do not enter AMI_Validator, but use a reasonable x bound
    % so the canonical config remains well formed.
    cfg = build_fso_config(M,P_avg,SNR_dB,sigma_X_sq, ...
        'ghN_h',40,'xMaxBound',max(5,2*P_avg),'saMaxIter',1, ...
        'saNStarts',2,'saUseParallel',false,'minGap',0,'pinZero',true, ...
        'seedInit',1,'logEvery',0);
end


function S = build_summary(caseGrid,MVec,gapVec,bitsVec,repVec,thresholds)
    nM = numel(MVec); nG = numel(gapVec); nB = numel(bitsVec); nR = numel(repVec);
    sz4 = [nM nG nB nR];
    sz3 = [nM nG nR];

    S.sourceAMI = nan(sz3);
    S.sourceRevalidationGap = nan(sz3);
    S.quantizedAMI = nan(sz4);
    S.amiPenalty = nan(sz4);
    S.deltaAMI = nan(sz4);
    S.absAMIChange = nan(sz4);
    S.quantizedShapingGain = nan(sz4);
    S.uniqueLevels = nan(sz4);
    S.collisionCount = nan(sz4);
    S.minAdjacentSpacing = nan(sz4);
    S.minDistinctSpacing = nan(sz4);
    S.quantStep = nan(sz4);
    S.xMax = nan(sz4);
    S.rmsDisplacement = nan(sz4);
    S.maxAbsDisplacement = nan(sz4);
    S.powerCalibrationScale = nan(sz4);
    S.meanPowerErrorBeforeCalibration = nan(sz4);
    S.meanPowerErrorAfterCalibration = nan(sz4);
    S.clippingCount = nan(sz4);

    for iM = 1:nM
        for iG = 1:nG
            for iR = 1:nR
                c = caseGrid{iM,iG,iR};
                S.sourceAMI(iM,iG,iR) = c.source.revalidatedAMI;
                S.sourceRevalidationGap(iM,iG,iR) = c.source.revalidationGap;
                for iB = 1:nB
                    q = c.quantized(iB);
                    S.quantizedAMI(iM,iG,iB,iR) = q.amiValidated;
                    S.amiPenalty(iM,iG,iB,iR) = q.amiPenalty;
                    S.deltaAMI(iM,iG,iB,iR) = q.deltaAMI;
                    S.absAMIChange(iM,iG,iB,iR) = q.absAMIChange;
                    S.quantizedShapingGain(iM,iG,iB,iR) = q.quantizedShapingGain;
                    S.uniqueLevels(iM,iG,iB,iR) = q.nUniqueLevels;
                    S.collisionCount(iM,iG,iB,iR) = q.collisionCount;
                    S.minAdjacentSpacing(iM,iG,iB,iR) = q.minAdjacentSpacing;
                    S.minDistinctSpacing(iM,iG,iB,iR) = q.minDistinctSpacing;
                    S.quantStep(iM,iG,iB,iR) = q.quantStepAfterCalibration;
                    S.xMax(iM,iG,iB,iR) = q.xMax;
                    S.rmsDisplacement(iM,iG,iB,iR) = q.rmsDisplacement;
                    S.maxAbsDisplacement(iM,iG,iB,iR) = q.maxAbsDisplacement;
                    S.powerCalibrationScale(iM,iG,iB,iR) = q.powerCalibrationScale;
                    S.meanPowerErrorBeforeCalibration(iM,iG,iB,iR) = q.meanPowerErrorBeforeCalibration;
                    S.meanPowerErrorAfterCalibration(iM,iG,iB,iR) = q.meanPowerErrorAfterCalibration;
                    S.clippingCount(iM,iG,iB,iR) = q.clippingCount;
                end
            end
        end
    end

    fields4 = {'quantizedAMI','amiPenalty','deltaAMI','absAMIChange','quantizedShapingGain', ...
        'uniqueLevels','collisionCount','minAdjacentSpacing','minDistinctSpacing','quantStep', ...
        'xMax','rmsDisplacement','maxAbsDisplacement','powerCalibrationScale', ...
        'meanPowerErrorBeforeCalibration','meanPowerErrorAfterCalibration','clippingCount'};
    S.mean = struct();
    S.stdAcrossReplicates = struct();
    for k = 1:numel(fields4)
        f = fields4{k};
        S.mean.(f) = mean(S.(f),4,'omitnan');
        S.stdAcrossReplicates.(f) = std(S.(f),0,4,'omitnan');
    end
    S.mean.sourceAMI = mean(S.sourceAMI,3,'omitnan');
    S.stdAcrossReplicates.sourceAMI = std(S.sourceAMI,0,3,'omitnan');

    S.allDistinctFraction = nan(nM,nG,nB);
    S.maxAbsAMIChangeAcrossReplicates = nan(nM,nG,nB);
    S.maxAMIPenaltyAcrossReplicates = nan(nM,nG,nB);
    for iM = 1:nM
        for iG = 1:nG
            for iB = 1:nB
                valsUnique = reshape(S.uniqueLevels(iM,iG,iB,:),1,[]);
                S.allDistinctFraction(iM,iG,iB) = mean(valsUnique == MVec(iM));
                valsAbs = reshape(S.absAMIChange(iM,iG,iB,:),1,[]);
                valsPenalty = reshape(S.amiPenalty(iM,iG,iB,:),1,[]);
                S.maxAbsAMIChangeAcrossReplicates(iM,iG,iB) = max(valsAbs,[],'omitnan');
                S.maxAMIPenaltyAcrossReplicates(iM,iG,iB) = max(valsPenalty,[],'omitnan');
            end
        end
    end

    S.requirements = struct();
    S.requirements.minBitsAllDistinct = nan(nM,nG);
    for iM = 1:nM
        for iG = 1:nG
            idx = find(reshape(S.allDistinctFraction(iM,iG,:),1,[]) == 1,1,'first');
            if ~isempty(idx), S.requirements.minBitsAllDistinct(iM,iG) = bitsVec(idx); end
        end
    end

    S.requirements.lossThresholds = thresholds;
    S.requirements.minBitsAbsAMIChange = nan(nM,nG,numel(thresholds));
    for it = 1:numel(thresholds)
        th = thresholds(it);
        for iM = 1:nM
            for iG = 1:nG
                vals = reshape(S.maxAbsAMIChangeAcrossReplicates(iM,iG,:),1,[]);
                idx = find(vals <= th,1,'first');
                if ~isempty(idx), S.requirements.minBitsAbsAMIChange(iM,iG,it) = bitsVec(idx); end
            end
        end
    end
end


function files = make_summary_plots(bundle,figDir,saveFigures)
    MVec = bundle.axes.M; gapVec = bundle.axes.minGap; bitsVec = bundle.axes.bits;
    S = bundle.summary; files = {};

    for iM = 1:numel(MVec)
        M = MVec(iM);
        f1 = figure('Name',sprintf('Quantization sensitivity M=%d',M),'Color','w');
        tl = tiledlayout(f1,1,2,'TileSpacing','compact','Padding','compact');

        ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
        for iG = 1:numel(gapVec)
            y = reshape(S.mean.amiPenalty(iM,iG,:),1,[]);
            e = reshape(S.stdAcrossReplicates.amiPenalty(iM,iG,:),1,[]);
            errorbar(ax1,bitsVec,y,e,'-o','LineWidth',1.3, ...
                'DisplayName',sprintf('d_{min}=%.4g',gapVec(iG)));
        end
        yline(ax1,0,'--');
        xlabel(ax1,'DAC resolution B [bits]');
        ylabel(ax1,'Validated AMI penalty [bits/symbol]');
        title(ax1,sprintf('Quantization penalty, M=%d',M));
        legend(ax1,'Location','best');

        ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
        for iG = 1:numel(gapVec)
            y = reshape(S.mean.uniqueLevels(iM,iG,:),1,[]);
            plot(ax2,bitsVec,y,'-o','LineWidth',1.3, ...
                'DisplayName',sprintf('d_{min}=%.4g',gapVec(iG)));
        end
        yline(ax2,M,'--','M distinct levels');
        xlabel(ax2,'DAC resolution B [bits]');
        ylabel(ax2,'Mean number of unique transmitted levels');
        title(ax2,'DAC-code collisions');
        legend(ax2,'Location','best');

        if saveFigures
            files = [files save_figure_pair(f1,figDir,sprintf('quantization_penalty_M%d',M))]; %#ok<AGROW>
        end

        f2 = figure('Name',sprintf('Quantized shaping gain M=%d',M),'Color','w');
        ax = axes(f2); hold(ax,'on'); grid(ax,'on');
        for iG = 1:numel(gapVec)
            y = reshape(S.mean.quantizedShapingGain(iM,iG,:),1,[]);
            e = reshape(S.stdAcrossReplicates.quantizedShapingGain(iM,iG,:),1,[]);
            errorbar(ax,bitsVec,y,e,'-o','LineWidth',1.3, ...
                'DisplayName',sprintf('d_{min}=%.4g',gapVec(iG)));
        end
        xlabel(ax,'DAC resolution B [bits]');
        ylabel(ax,'Validated quantized GS - quantized PAM [bits/symbol]');
        title(ax,sprintf('Quantized shaping gain: M=%d, SNR=%.1f dB, \\sigma_X^2=%.3f', ...
            M,bundle.axes.SNR_dB,bundle.axes.sigma_X_sq));
        legend(ax,'Location','best');
        if saveFigures
            files = [files save_figure_pair(f2,figDir,sprintf('quantized_gain_M%d',M))]; %#ok<AGROW>
        end
    end
end


function print_summary(bundle)
    MVec = bundle.axes.M; gapVec = bundle.axes.minGap; bitsVec = bundle.axes.bits;
    S = bundle.summary;
    fprintf('\nQuantization sensitivity summary\n');
    fprintf('     M |    d_min | source AMI | min bits distinct | min bits |dAMI|<=1e-2 | min bits |dAMI|<=1e-3\n');
    fprintf('------------------------------------------------------------------------------------------------\n');
    for iM = 1:numel(MVec)
        for iG = 1:numel(gapVec)
            bDistinct = S.requirements.minBitsAllDistinct(iM,iG);
            b1 = threshold_bits(S,1e-2,iM,iG);
            b2 = threshold_bits(S,1e-3,iM,iG);
            fprintf('%6d | %8.4f | %10.6f | %17s | %21s | %21s\n', ...
                MVec(iM),gapVec(iG),S.mean.sourceAMI(iM,iG), ...
                bit_text(bDistinct),bit_text(b1),bit_text(b2));
        end
    end

    fprintf('\nMean AMI penalty at each DAC resolution:\n');
    fprintf('bits = %s\n',mat2str(bitsVec));
    for iM = 1:numel(MVec)
        for iG = 1:numel(gapVec)
            vals = reshape(S.mean.amiPenalty(iM,iG,:),1,[]);
            fprintf('M=%d d_min=%.4g : %s\n',MVec(iM),gapVec(iG),mat2str(vals,5));
        end
    end

    fprintf('\nMax source-AMI revalidation mismatch: %.3e bits/symbol\n', ...
        max(abs(S.sourceRevalidationGap(:)),[],'omitnan'));
end


function b = threshold_bits(S,threshold,iM,iG)
    vals = S.requirements.lossThresholds;
    [dist,idx] = min(abs(vals-threshold));
    if isempty(idx) || dist > max(1e-12,1e-9*threshold)
        b = NaN;
    else
        b = S.requirements.minBitsAbsAMIChange(iM,iG,idx);
    end
end


function t = bit_text(v)
    if isnan(v), t = 'not reached'; else, t = sprintf('%g',v); end
end


function [src,description] = load_source_bundle(source)
    if isstruct(source)
        src = source;
        description = '<workspace bundle>';
        return;
    end

    path = char(string(source));
    if ~isfile(path)
        error('sim_quantization_sensitivity:SourceNotFound','Source file not found: %s',path);
    end
    data = load(path);
    if isfield(data,'bundle')
        src = data.bundle;
    elseif isfield(data,'qbundle')
        error('sim_quantization_sensitivity:WrongSourceType', ...
            'The supplied file contains a quantization bundle, not a minimum-spacing source bundle.');
    else
        names = fieldnames(data);
        if numel(names)==1 && isstruct(data.(names{1}))
            src = data.(names{1});
        else
            error('sim_quantization_sensitivity:BundleMissing', ...
                'Could not identify a source bundle in %s.',path);
        end
    end
    description = path;
end


function validate_source_bundle(src)
    requiredTop = {'axes','summary'};
    for k = 1:numel(requiredTop)
        if ~isfield(src,requiredTop{k})
            error('sim_quantization_sensitivity:BadSource', ...
                'Source bundle is missing field %s.',requiredTop{k});
        end
    end
    requiredAxes = {'M','minGap','replicate','SNR_dB','sigma_X_sq','P_avg'};
    for k = 1:numel(requiredAxes)
        if ~isfield(src.axes,requiredAxes{k})
            error('sim_quantization_sensitivity:BadSource', ...
                'Source axes missing field %s.',requiredAxes{k});
        end
    end
    if ~isfield(src.summary,'xWinner') || ~isfield(src.summary,'amiGSValidated')
        error('sim_quantization_sensitivity:BadSource', ...
            'Source summary must contain xWinner and amiGSValidated.');
    end
end


function idx = map_axis_values(requested,available,name)
    idx = zeros(size(requested));
    scale = max(1,max(abs(available)));
    tol = 1e-10*scale;
    for k = 1:numel(requested)
        j = find(abs(available-requested(k)) <= tol,1,'first');
        if isempty(j)
            error('sim_quantization_sensitivity:AxisValueMissing', ...
                'Requested %s value %.12g is not present in source bundle.',name,requested(k));
        end
        idx(k) = j;
    end
end


function id = source_run_id(src)
    id = '';
    if isfield(src,'run') && isfield(src.run,'id'), id = src.run.id; end
end


function name = source_experiment(src)
    name = '';
    if isfield(src,'experiment'), name = src.experiment; end
end


function d = safe_min_diff(x)
    x = sort(real(x(:)),'ascend');
    if numel(x) < 2, d = NaN; else, d = min(diff(x)); end
end


function files = save_figure_pair(fig,figDir,baseName)
    pngPath = fullfile(figDir,[baseName '.png']);
    figPath = fullfile(figDir,[baseName '.fig']);
    exportgraphics(fig,pngPath,'Resolution',200);
    savefig(fig,figPath);
    files = {pngPath,figPath};
end


function token = sanitize_run_label(v)
    token = strtrim(char(string(v)));
    if isempty(token), token = 'run'; return; end
    token = regexprep(token,'[^A-Za-z0-9._-]+','_');
    token = regexprep(token,'_+','_');
    token = regexprep(token,'^[._-]+|[._-]+$','');
    if isempty(token), token = 'run'; end
end


function [outDir,runId] = make_unique_run_directory(conditionDir,requestedRunId)
    if ~exist(conditionDir,'dir'), mkdir(conditionDir); end
    runId = requestedRunId;
    outDir = fullfile(conditionDir,runId);
    suffix = 1;
    while exist(outDir,'dir')
        suffix = suffix+1;
        runId = sprintf('%s_%02d',requestedRunId,suffix);
        outDir = fullfile(conditionDir,runId);
    end
    [ok,msg] = mkdir(outDir);
    if ~ok
        error('sim_quantization_sensitivity:RunDirectoryCreateFailed', ...
            'Could not create run directory %s: %s',outDir,msg);
    end
end


function token = join_numeric_token(v)
    parts = arrayfun(@(x) sprintf('%d',x),v,'UniformOutput',false);
    token = strjoin(parts,'-');
end


function token = num_token(v,nDigits)
    if nargin < 2, nDigits = 4; end
    fmt = sprintf('%%+.%df',nDigits);
    token = sprintf(fmt,v);
    token = strrep(token,'+','p');
    token = strrep(token,'-','m');
    token = strrep(token,'.','p');
end
