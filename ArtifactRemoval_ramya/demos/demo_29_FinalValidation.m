clear; close all; clc;

%% demo_29_FinalValidation_WithSMARTAFull
%
% Purpose:
% Final real-data validation across the selected 31 V1 electrodes,
% now including true optimal-shrinkage SMARTAFull.
%
% Methods:
%   1. Raw high-stimulation
%   2. ERPSubtraction
%   3. PCATemplate K=10
%   4. SMARTALite ensemble K3K5
%   5. SMARTAFull K=3, hp=100 Hz
%
% Main outputs:
%   demo29_final_validation_long.csv
%   demo29_final_validation_summary.csv
%   demo29_best_method_counts.csv
%
% Important:
% The balanced score here is a normalized per-electrode score across
% candidate cleaning methods. Lower is better. Do not directly mix this
% balanced-score value with older scripts if the older formula was different.

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName  = 'dona';
gridType     = 'Microelectrode';
expDate      = '290825';
protocolName = 'GRF_001';

% Conditions used in previous final validation
noStimCondition   = {1,1,1,5,5,4};
highStimCondition = {7,1,1,5,5,4};

% 31 selected V1 electrodes used in final validation
goodV1Electrodes = [ ...
     3  4  5  6  7  8 10 11 13 ...
    14 15 17 18 19 20 21 22 23 ...
    24 25 26 27 28 29 30 31 34 ...
    35 36 41 42];

% Analysis windows
artifactWindow = [0 0.4];
fftWindow      = [0 0.4];

% Frequency metrics
freqRangeForMetric = [0 200];
stimFreq = 20;
harmonicsToUse = stimFreq:stimFreq:200;
harmonicHalfBandwidthHz = 1.0;

% ICMS pulse settings
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Methods
methodList = { ...
    'RawHighStim', ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_Ensemble_K3K5', ...
    'SMARTAFull_K3_hp100'};

candidateMethods = { ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_Ensemble_K3K5', ...
    'SMARTAFull_K3_hp100'};

saveFigures = true;

%% ========================================================================
% Paths
% ========================================================================

projectMethods = fullfile(folderSourceString,'ICMS_Artifact_Removal','methods');
projectUtils   = fullfile(folderSourceString,'ICMS_Artifact_Removal','utils');
githubSMARTA   = fullfile(folderSourceString,'artifact_removal','Artifact-removal-main');

% GitHub folder supplies createPseudoNoise.m.
% Project methods folder should remain first for SMARTAFull and
% optimal_shrinkage_color_fast.
if exist(githubSMARTA,'dir')
    addpath(githubSMARTA);
end

if exist(projectUtils,'dir')
    addpath(projectUtils);
end

if exist(projectMethods,'dir')
    addpath(projectMethods,'-begin');
end

fprintf('Checking active SMARTA paths:\n');
disp(which('SMARTAFull','-all'));
disp(which('optimal_shrinkage_color_fast','-all'));
disp(which('createPseudoNoise','-all'));

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder   = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

resultFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','metrics','demo29_FinalValidation_WithSMARTAFull');

figFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','figures','demo29_FinalValidation_WithSMARTAFull');

if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Load shared files
% ========================================================================

I = load(lfpInfoFile);
P = load(paramFile);

timeVals = I.timeVals(:).';
parameterCombinations = P.parameterCombinations;

noStimTrials   = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

fprintf('\nNo-stim trials   : %d\n',numel(noStimTrials));
fprintf('High-stim trials : %d\n',numel(highStimTrials));
fprintf('Electrodes       : %d\n',numel(goodV1Electrodes));

fftIdx = find(timeVals > fftWindow(1) & timeVals <= fftWindow(2));

Fs = 1/median(diff(timeVals));
nFFT = length(fftIdx);
freqAxis = (0:nFFT-1) * (Fs/nFFT);

freqMask = freqAxis >= freqRangeForMetric(1) & freqAxis <= freqRangeForMetric(2);
harmonicMask = buildHarmonicMask(freqAxis,harmonicsToUse,harmonicHalfBandwidthHz,freqRangeForMetric);

fprintf('FFT samples      : %d\n',nFFT);
fprintf('Frequency step   : %.2f Hz\n',Fs/nFFT);

%% ========================================================================
% Method options
% ========================================================================

opts = struct();
opts.artifactWindow = artifactWindow;
opts.pulseTimes = pulseTimes;
opts.stiFreq = smartaStimFreq;
opts.prePulse = smartaPrePulse;

%% ========================================================================
% Main validation loop
% ========================================================================

rows = {};
rowCounter = 0;

for iElec = 1:numel(goodV1Electrodes)

    elecNum = goodV1Electrodes(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    fprintf('\n============================================================\n');
    fprintf('Processing elec%d (%d/%d)\n',elecNum,iElec,numel(goodV1Electrodes));
    fprintf('============================================================\n');

    if ~exist(elecFile,'file')
        warning('Missing %s. Skipping.',elecFile);
        continue;
    end

    D = load(elecFile);

    if ~isfield(D,'analogData')
        warning('analogData missing in elec%d. Skipping.',elecNum);
        continue;
    end

    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    fftNoStim = computeLogFFT(dataNoStim,fftIdx);
    trialSTDNoStim = std(dataNoStim(:,fftIdx),0,1);
    noStimRMS = safeRMS(dataNoStim(:));

    for iMethod = 1:numel(methodList)

        methodName = methodList{iMethod};

        fprintf('Method: %s\n',methodName);

        if strcmp(methodName,'RawHighStim')
            cleanedHigh = dataHighStim;
            cleanedNoStim = dataNoStim;
            usedShrinkageHigh = NaN;
            shrinkageErrorHigh = '';
            usedShrinkageNoStim = NaN;
            shrinkageErrorNoStim = '';
        else
            [cleanedHigh,diagHigh] = applyCleaningMethod(dataHighStim,timeVals,methodName,opts);
            [cleanedNoStim,diagNoStim] = applyCleaningMethod(dataNoStim,timeVals,methodName,opts);

            usedShrinkageHigh = diagHigh.usedShrinkage;
            shrinkageErrorHigh = diagHigh.shrinkageErrorMessage;
            usedShrinkageNoStim = diagNoStim.usedShrinkage;
            shrinkageErrorNoStim = diagNoStim.shrinkageErrorMessage;
        end

        fftCleanHigh = computeLogFFT(cleanedHigh,fftIdx);
        fftCleanNoStim = computeLogFFT(cleanedNoStim,fftIdx);

        diffHigh = fftCleanHigh - fftNoStim;

        totalFFTError = norm(diffHigh(freqMask));
        harmonicExcessError = norm(max(diffHigh(harmonicMask),0));
        harmonicAbsError = norm(diffHigh(harmonicMask));
        belowReferenceError = norm(max(-diffHigh(freqMask),0));

        trialSTDCleanHigh = std(cleanedHigh(:,fftIdx),0,1);
        trialSTDError = safeRMS(trialSTDCleanHigh - trialSTDNoStim);
        trialSTDNormError = trialSTDError / max(safeRMS(trialSTDNoStim),eps);

        noStimDiff = cleanedNoStim - dataNoStim;
        noStimRMSDistortion = safeRMS(noStimDiff(:));
        noStimNormRMSDistortion = noStimRMSDistortion / max(noStimRMS,eps);
        noStimFFTDistortion = norm(fftCleanNoStim(freqMask) - fftNoStim(freqMask));

        highMeanTrace = mean(cleanedHigh,1);
        refMeanTrace = mean(dataNoStim,1);
        meanTraceCorrelation = safeCorr(highMeanTrace(fftIdx).',refMeanTrace(fftIdx).');

        rowCounter = rowCounter + 1;

        rows(rowCounter,:) = { ...
            elecNum, ...
            methodName, ...
            totalFFTError, ...
            harmonicExcessError, ...
            harmonicAbsError, ...
            belowReferenceError, ...
            trialSTDError, ...
            trialSTDNormError, ...
            noStimRMSDistortion, ...
            noStimNormRMSDistortion, ...
            noStimFFTDistortion, ...
            meanTraceCorrelation, ...
            usedShrinkageHigh, ...
            string(shrinkageErrorHigh), ...
            usedShrinkageNoStim, ...
            string(shrinkageErrorNoStim)};
    end
end

T = cell2table(rows, ...
    'VariableNames',{ ...
    'electrode', ...
    'method', ...
    'totalFFTError', ...
    'harmonicExcessError', ...
    'harmonicAbsError', ...
    'belowReferenceError', ...
    'trialSTDError', ...
    'trialSTDNormError', ...
    'noStimRMSDistortion', ...
    'noStimNormRMSDistortion', ...
    'noStimFFTDistortion', ...
    'meanTraceCorrelation', ...
    'usedShrinkageHigh', ...
    'shrinkageErrorHigh', ...
    'usedShrinkageNoStim', ...
    'shrinkageErrorNoStim'});

%% ========================================================================
% Normalized balanced score
% ========================================================================

T.balancedScore = nan(height(T),1);

balancedComponents = { ...
    'totalFFTError', ...
    'harmonicExcessError', ...
    'belowReferenceError', ...
    'trialSTDNormError', ...
    'noStimNormRMSDistortion', ...
    'noStimFFTDistortion'};

uniqueElectrodes = unique(T.electrode);

for iElec = 1:numel(uniqueElectrodes)

    elecNum = uniqueElectrodes(iElec);

    idxCandidate = T.electrode == elecNum & ismember(T.method,candidateMethods);

    if ~any(idxCandidate)
        continue;
    end

    score = zeros(sum(idxCandidate),1);

    for iComp = 1:numel(balancedComponents)

        compName = balancedComponents{iComp};
        vals = T.(compName)(idxCandidate);

        minVal = min(vals,[],'omitnan');
        maxVal = max(vals,[],'omitnan');

        if ~isfinite(minVal) || ~isfinite(maxVal) || abs(maxVal-minVal) < eps
            normVals = zeros(size(vals));
        else
            normVals = (vals - minVal) ./ (maxVal - minVal);
        end

        score = score + normVals;
    end

    score = score ./ numel(balancedComponents);

    tmpIdx = find(idxCandidate);
    T.balancedScore(tmpIdx) = score;
end

%% ========================================================================
% Summary tables
% ========================================================================

summaryTable = summarizeFinalTable(T,methodList);

disp(' ');
disp('============================================================');
disp('Final validation summary including SMARTAFull');
disp('============================================================');
disp(summaryTable);

writetable(T,fullfile(resultFolder,'demo29_final_validation_long.csv'));
writetable(summaryTable,fullfile(resultFolder,'demo29_final_validation_summary.csv'));

%% Best-method counts

bestCountsTable = computeBestMethodCounts(T,candidateMethods);

disp(' ');
disp('============================================================');
disp('Best-method counts');
disp('============================================================');
disp(bestCountsTable);

writetable(bestCountsTable,fullfile(resultFolder,'demo29_best_method_counts.csv'));

%% ========================================================================
% Save figures
% ========================================================================

if saveFigures

    fig1 = figure('Color','w','Position',[100 100 1050 500]);
    bar(summaryTable.meanTotalFFTError);
    set(gca,'XTick',1:height(summaryTable),'XTickLabel',prettyMethodNames(summaryTable.method));
    xtickangle(25);
    ylabel('Mean total FFT error');
    title('Final validation: total FFT error');
    grid on; box off;
    saveAllFigureFormats(fig1,figFolder,'demo29_mean_total_fft_error');

    fig2 = figure('Color','w','Position',[100 100 1050 500]);
    bar(summaryTable.meanHarmonicExcessError);
    set(gca,'XTick',1:height(summaryTable),'XTickLabel',prettyMethodNames(summaryTable.method));
    xtickangle(25);
    ylabel('Mean harmonic excess error');
    title('Final validation: stimulation-harmonic excess');
    grid on; box off;
    saveAllFigureFormats(fig2,figFolder,'demo29_mean_harmonic_excess_error');

    fig3 = figure('Color','w','Position',[100 100 1050 500]);
    bar(summaryTable.meanNoStimFFTDistortion);
    set(gca,'XTick',1:height(summaryTable),'XTickLabel',prettyMethodNames(summaryTable.method));
    xtickangle(25);
    ylabel('Mean no-stimulation FFT distortion');
    title('Negative-control distortion');
    grid on; box off;
    saveAllFigureFormats(fig3,figFolder,'demo29_mean_no_stim_fft_distortion');

    fig4 = figure('Color','w','Position',[100 100 1050 500]);
    bar(summaryTable.meanBalancedScore);
    set(gca,'XTick',1:height(summaryTable),'XTickLabel',prettyMethodNames(summaryTable.method));
    xtickangle(25);
    ylabel('Mean normalized balanced score');
    title('Final validation: normalized balanced score');
    grid on; box off;
    saveAllFigureFormats(fig4,figFolder,'demo29_mean_balanced_score');

    fig5 = figure('Color','w','Position',[100 100 1050 500]);
    bar(summaryTable.meanTrialSTDNormError);
    set(gca,'XTick',1:height(summaryTable),'XTickLabel',prettyMethodNames(summaryTable.method));
    xtickangle(25);
    ylabel('Mean normalized trial-STD error');
    title('Trial-variability preservation');
    grid on; box off;
    saveAllFigureFormats(fig5,figFolder,'demo29_mean_trial_std_norm_error');
end

%% ========================================================================
% Final printed interpretation helper
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_29 complete\n');
fprintf('============================================================\n');
fprintf('Saved tables in:\n%s\n',resultFolder);
fprintf('\nSaved figures in:\n%s\n',figFolder);

fprintf('\nKey files:\n');
fprintf('  demo29_final_validation_long.csv\n');
fprintf('  demo29_final_validation_summary.csv\n');
fprintf('  demo29_best_method_counts.csv\n');

fprintf('\nNext: send demo29_final_validation_summary.csv and best counts.\n');

%% ========================================================================
% Helper functions
% ========================================================================

function [signalClean,diagOut] = applyCleaningMethod(signalRaw,timeVals,methodName,opts)

    diagOut = struct();
    diagOut.usedShrinkage = NaN;
    diagOut.shrinkageErrorMessage = "";

    switch methodName

        case 'ERPSubtraction'

            erpParams = struct();
            erpParams.subtractWindow = opts.artifactWindow;
            erpParams.doBaselineCorrection = false;

            out = ERPSubtraction(signalRaw,timeVals,erpParams);
            signalClean = getCleanedData(out);

        case 'PCATemplate_K10'

            pcaParams = struct();
            pcaParams.artifactWindow = opts.artifactWindow;
            pcaParams.numComponents = 10;
            pcaParams.removeMeanTemplate = true;
            pcaParams.taperEdgeMS = 2;
            pcaParams.doBaselineCorrection = false;

            out = PCATemplate(signalRaw,timeVals,pcaParams);
            signalClean = getCleanedData(out);

        case 'SMARTALite_Ensemble_K3K5'

            liteOpts = struct();
            liteOpts.pulseTimes = opts.pulseTimes;
            liteOpts.stiFreq = opts.stiFreq;
            liteOpts.prePulse = opts.prePulse;
            liteOpts.window = opts.artifactWindow;
            liteOpts.computeFFT = false;
            liteOpts.excludeSameTrial = true;

            liteOpts.K = 3;
            outK3 = SMARTALite(signalRaw,timeVals,liteOpts);
            cleanK3 = getCleanedData(outK3);

            liteOpts.K = 5;
            outK5 = SMARTALite(signalRaw,timeVals,liteOpts);
            cleanK5 = getCleanedData(outK5);

            signalClean = 0.5 * (cleanK3 + cleanK5);

        case 'SMARTAFull_K3_hp100'

            fullOpts = struct();
            fullOpts.pulseTimes = opts.pulseTimes;
            fullOpts.stiFreq = opts.stiFreq;
            fullOpts.prePulse = opts.prePulse;
            fullOpts.window = opts.artifactWindow;
            fullOpts.computeFFT = false;
            fullOpts.excludeSameTrial = true;

            fullOpts.K = 3;
            fullOpts.hpCutoff = 100;
            fullOpts.shrinkLoss = 'fro';
            fullOpts.shrinkKL = 10;
            fullOpts.shrinkKH = 15;

            out = SMARTAFull(signalRaw,timeVals,fullOpts);
            signalClean = getCleanedData(out);

            if isfield(out,'diagnostics')
                if isfield(out.diagnostics,'usedShrinkage')
                    diagOut.usedShrinkage = out.diagnostics.usedShrinkage;
                end

                if isfield(out.diagnostics,'shrinkageErrorMessage')
                    diagOut.shrinkageErrorMessage = string(out.diagnostics.shrinkageErrorMessage);
                end
            end

        otherwise
            error('Unknown methodName: %s',methodName);
    end
end

function cleanedData = getCleanedData(methodOut)

    if isfield(methodOut,'cleanedData')
        cleanedData = methodOut.cleanedData;
    elseif isfield(methodOut,'cleanedTrials')
        cleanedData = methodOut.cleanedTrials;
    elseif isfield(methodOut,'z')
        cleanedData = methodOut.z;
    else
        error('Method output has no cleaned data field.');
    end
end

function y = computeLogFFT(data,fftIdx)

    fftVals = abs(fft(data(:,fftIdx),[],2));
    y = log10(mean(fftVals,1) + eps);
end

function harmonicMask = buildHarmonicMask(freqAxis,harmonicsToUse,halfBandwidthHz,freqRangeForMetric)

    harmonicMask = false(size(freqAxis));

    for iH = 1:numel(harmonicsToUse)

        h = harmonicsToUse(iH);

        if h < freqRangeForMetric(1) || h > freqRangeForMetric(2)
            continue;
        end

        harmonicMask = harmonicMask | abs(freqAxis - h) <= halfBandwidthHz;
    end
end

function r = safeRMS(x)

    x = x(:);
    r = sqrt(mean(x.^2,'omitnan'));
end

function c = safeCorr(x,y)

    x = x(:);
    y = y(:);

    good = isfinite(x) & isfinite(y);

    if sum(good) < 3
        c = NaN;
        return;
    end

    x = x(good);
    y = y(good);

    if std(x) == 0 || std(y) == 0
        c = NaN;
        return;
    end

    C = corrcoef(x,y);
    c = C(1,2);
end

function summaryTable = summarizeFinalTable(T,methodList)

    rows = {};

    for iMethod = 1:numel(methodList)

        methodName = methodList{iMethod};
        idx = strcmp(T.method,methodName);

        rows(end+1,:) = { ...
            methodName, ...
            sum(idx), ...
            mean(T.totalFFTError(idx),'omitnan'), ...
            std(T.totalFFTError(idx),'omitnan'), ...
            mean(T.harmonicExcessError(idx),'omitnan'), ...
            std(T.harmonicExcessError(idx),'omitnan'), ...
            mean(T.harmonicAbsError(idx),'omitnan'), ...
            mean(T.belowReferenceError(idx),'omitnan'), ...
            mean(T.trialSTDNormError(idx),'omitnan'), ...
            mean(T.noStimNormRMSDistortion(idx),'omitnan'), ...
            mean(T.noStimFFTDistortion(idx),'omitnan'), ...
            mean(T.meanTraceCorrelation(idx),'omitnan'), ...
            mean(T.balancedScore(idx),'omitnan'), ...
            std(T.balancedScore(idx),'omitnan'), ...
            mean(T.usedShrinkageHigh(idx),'omitnan')};
    end

    summaryTable = cell2table(rows, ...
        'VariableNames',{ ...
        'method', ...
        'nElectrodes', ...
        'meanTotalFFTError', ...
        'stdTotalFFTError', ...
        'meanHarmonicExcessError', ...
        'stdHarmonicExcessError', ...
        'meanHarmonicAbsError', ...
        'meanBelowReferenceError', ...
        'meanTrialSTDNormError', ...
        'meanNoStimNormRMSDistortion', ...
        'meanNoStimFFTDistortion', ...
        'meanTraceCorrelation', ...
        'meanBalancedScore', ...
        'stdBalancedScore', ...
        'meanUsedShrinkageHigh'});
end

function bestCountsTable = computeBestMethodCounts(T,candidateMethods)

    metrics = { ...
        'totalFFTError', ...
        'harmonicExcessError', ...
        'harmonicAbsError', ...
        'belowReferenceError', ...
        'trialSTDNormError', ...
        'noStimNormRMSDistortion', ...
        'noStimFFTDistortion', ...
        'balancedScore'};

    rows = {};
    uniqueElectrodes = unique(T.electrode);

    for iMetric = 1:numel(metrics)

        metricName = metrics{iMetric};

        counts = zeros(numel(candidateMethods),1);

        for iElec = 1:numel(uniqueElectrodes)

            elecNum = uniqueElectrodes(iElec);
            idx = T.electrode == elecNum & ismember(T.method,candidateMethods);

            if ~any(idx)
                continue;
            end

            vals = T.(metricName)(idx);
            methodsHere = T.method(idx);

            [~,bestIdx] = min(vals);

            bestMethod = methodsHere{bestIdx};
            countIdx = find(strcmp(candidateMethods,bestMethod));

            if ~isempty(countIdx)
                counts(countIdx) = counts(countIdx) + 1;
            end
        end

        for iMethod = 1:numel(candidateMethods)

            rows(end+1,:) = { ...
                metricName, ...
                candidateMethods{iMethod}, ...
                counts(iMethod)};
        end
    end

    bestCountsTable = cell2table(rows, ...
        'VariableNames',{'metric','method','bestCount'});
end

function names = prettyMethodNames(methodList)

    names = cell(size(methodList));

    for i = 1:numel(methodList)
        names{i} = prettyMethodName(methodList{i});
    end
end

function name = prettyMethodName(methodName)

    if isstring(methodName)
        methodName = char(methodName);
    end

    switch methodName
        case 'RawHighStim'
            name = 'Raw high-stim';
        case 'ERPSubtraction'
            name = 'ERPSubtraction';
        case 'PCATemplate_K10'
            name = 'PCATemplate K=10';
        case 'SMARTALite_Ensemble_K3K5'
            name = 'SMARTALite ensemble K3K5';
        case 'SMARTAFull_K3_hp100'
            name = 'SMARTAFull K=3, hp=100 Hz';
        otherwise
            name = methodName;
    end
end

function saveAllFigureFormats(figHandle,figFolder,baseName)

    pngFile = fullfile(figFolder,[baseName '.png']);
    figFile = fullfile(figFolder,[baseName '.fig']);
    pdfFile = fullfile(figFolder,[baseName '.pdf']);

    saveas(figHandle,pngFile);
    savefig(figHandle,figFile);

    try
        exportgraphics(figHandle,pdfFile,'ContentType','vector');
    catch
        exportgraphics(figHandle,pdfFile,'ContentType','image');
    end

    fprintf('Saved figure:\n%s\n%s\n%s\n',pngFile,figFile,pdfFile);
end