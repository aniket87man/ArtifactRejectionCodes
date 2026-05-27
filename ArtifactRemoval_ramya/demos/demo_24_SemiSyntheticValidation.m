clear; close all; clc;

%% demo_24_SemiSyntheticValidation
%
% Goal:
% Validate artifact-removal methods using semi-synthetic data.
%
% Why:
% In real high-stim data, the true clean LFP is unknown.
% Here, we create semi-synthetic contaminated data:
%
%   syntheticData = cleanNoStimData + estimatedArtifactFromHighStim
%
% Therefore, the true clean target is known:
%
%   cleanTarget = cleanNoStimData
%
% Methods compared:
%   1. Raw semi-synthetic contaminated data
%   2. ERPSubtraction
%   3. PCATemplate K=20
%   4. SMARTALite K=3
%   5. SMARTALite K=5
%   6. SMARTALite Ensemble K3K5
%
% Main metrics:
%   1. Normalized time-domain recovery error
%   2. FFT recovery error vs true clean no-stim
%   3. Above-reference spectral error
%   4. Below-reference spectral error
%   5. Balanced spectral score
%   6. Harmonic excess error
%   7. RMS recovery error
%
% Lower is better for all error metrics.
%
% Important:
% This is not meant to replace real-data validation.
% It is a control analysis showing whether the methods can recover
% a known clean signal when artifact-like waveforms are added.

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName  = 'dona';
gridType     = 'Microelectrode';
expDate      = '290825';
protocolName = 'GRF_001';

% Conditions
noStimCondition   = {1,1,1,5,5,4};
highStimCondition = {7,1,1,5,5,4};

% Electrode selection
v1Electrodes = 1:48;
stimElectrode = 1;
excludeStimElectrode = true;

maxElectrodesToRun = Inf;

% Windows
artifactWindow = [-0.02 0.4];
fftWindow = [0 0.4];
timeMetricWindow = [0 0.4];
freqRangeForMetric = [0 200];

% Stimulation harmonic settings
stimFreq = 20;
harmonicsToUse = stimFreq:stimFreq:200;
harmonicHalfBandwidthHz = 1.0;

% SMARTALite full-cycle settings
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Semi-synthetic artifact settings
rngSeed = 24;
artifactScale = 1.0;

% Over-cleaning penalty
lambdaBelow = 1.5;

saveFigures = true;

%% ========================================================================
% Paths and loading
% ========================================================================

rng(rngSeed);

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder   = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

I = load(lfpInfoFile);
P = load(paramFile);

timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

noStimTrials = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));
timeMetricIdx = find(timeVals >= timeMetricWindow(1) & timeVals < timeMetricWindow(2));

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));
fprintf('FFT window samples = %d\n',length(fftIdx));
fprintf('FFT window duration = %.4f s\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% ========================================================================
% Check required functions
% ========================================================================

requiredMethods = { ...
    'ERPSubtraction', ...
    'PCATemplate', ...
    'SMARTALite', ...
    'compute_fft_summary'};

fprintf('\nChecking required functions:\n');

for i = 1:length(requiredMethods)
    fn = requiredMethods{i};
    fnPath = which(fn);

    if isempty(fnPath)
        warning('%s not found on MATLAB path.',fn);
    else
        fprintf('%s -> %s\n',fn,fnPath);
    end
end

%% ========================================================================
% Good V1 electrode selection
% ========================================================================

badChannels = [];

impedanceFileName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat');

if exist(impedanceFileName,'file')
    Z = load(impedanceFileName);

    if isfield(Z,'impedanceValues')
        impedanceValues = Z.impedanceValues;
    elseif isfield(Z,'electrodeImpedances')
        impedanceValues = Z.electrodeImpedances;
    else
        error('Could not find impedanceValues or electrodeImpedances in impedance file.');
    end

    badImpedanceCutoff = 2500;
    badChannels = unique([find(impedanceValues > badImpedanceCutoff), find(isnan(impedanceValues))]);
else
    warning('Could not find impedanceValues.mat. Using no impedance-based bad channels.');
end

rfDataFileName = [subjectName gridType 'RFData.mat'];
rfDataPath = which(rfDataFileName);

if isempty(rfDataPath)
    error('Could not find %s on MATLAB path.',rfDataFileName);
end

rfData = load(rfDataPath);

if ~isfield(rfData,'highRMSElectrodes')
    error('RFData file does not contain highRMSElectrodes.');
end

highRMSElectrodes = rfData.highRMSElectrodes(:)';

goodElectrodesAll = setdiff(highRMSElectrodes,badChannels);
goodV1Electrodes = intersect(goodElectrodesAll,v1Electrodes);

fprintf('\nGood V1 electrodes before stim-electrode exclusion:\n');
disp(goodV1Electrodes);
fprintf('Number of good V1 electrodes = %d\n',length(goodV1Electrodes));

if excludeStimElectrode
    goodV1Electrodes = setdiff(goodV1Electrodes,stimElectrode);
    fprintf('\nExcluding stimulation electrode elec%d.\n',stimElectrode);
end

fprintf('\nFinal good V1 electrodes used in demo_24:\n');
disp(goodV1Electrodes);
fprintf('Final number of electrodes = %d\n',length(goodV1Electrodes));

%% ========================================================================
% Electrode files
% ========================================================================

elecFiles = dir(fullfile(lfpFolder,'elec*.mat'));

elecNums = nan(length(elecFiles),1);

for i = 1:length(elecFiles)
    token = regexp(elecFiles(i).name,'^elec(\d+)\.mat$','tokens');

    if ~isempty(token)
        elecNums(i) = str2double(token{1}{1});
    end
end

validFileIdx = ~isnan(elecNums);
elecFiles = elecFiles(validFileIdx);
elecNums = elecNums(validFileIdx);

[elecNums,sortIdx] = sort(elecNums);
elecFiles = elecFiles(sortIdx);

keepIdx = ismember(elecNums,goodV1Electrodes);

elecNums = elecNums(keepIdx);
elecFiles = elecFiles(keepIdx);

if isfinite(maxElectrodesToRun)
    nToRun = min(maxElectrodesToRun,length(elecNums));
    elecNums = elecNums(1:nToRun);
    elecFiles = elecFiles(1:nToRun);
end

fprintf('\nNumber of electrode files selected for analysis: %d\n',length(elecNums));
fprintf('Electrode files selected:\n');
disp(elecNums');

%% ========================================================================
% Method parameters
% ========================================================================

% ERPSubtraction
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

% PCATemplate K=10
pcaParamsK10 = struct();
pcaParamsK10.artifactWindow = artifactWindow;
pcaParamsK10.numComponents = 10;
pcaParamsK10.removeMeanTemplate = true;
pcaParamsK10.taperEdgeMS = 2;
pcaParamsK10.doBaselineCorrection = false;

% SMARTALite K=3
smartaK3Opts = struct();
smartaK3Opts.pulseTimes = pulseTimes;
smartaK3Opts.stiFreq = smartaStimFreq;
smartaK3Opts.prePulse = smartaPrePulse;
smartaK3Opts.K = 3;
smartaK3Opts.window = fftWindow;
smartaK3Opts.computeFFT = false;

% SMARTALite K=5
smartaK5Opts = smartaK3Opts;
smartaK5Opts.K = 5;

methodNames = { ...
    'RawSemiSynthetic', ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_Ensemble_K3K5'};

candidateMethods = methodNames(2:end);

%% ========================================================================
% Output folders
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo24_SemiSyntheticValidation');

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Main loop
% ========================================================================

longRows = struct([]);
rowCounter = 0;

for iElec = 1:length(elecFiles)

    elecNum = elecNums(iElec);
    elecFile = fullfile(lfpFolder,elecFiles(iElec).name);

    fprintf('\n============================================================\n');
    fprintf('Processing electrode %d/%d: elec%d\n',iElec,length(elecFiles),elecNum);
    fprintf('============================================================\n');

    try
        D = load(elecFile);

        if ~isfield(D,'analogData')
            warning('Skipping elec%d: analogData not found.',elecNum);
            continue;
        end

        analogData = D.analogData;

        dataNoStimAll = analogData(noStimTrials,:);
        dataHighStimAll = analogData(highStimTrials,:);

        %% ----------------------------------------------------------------
        % Build semi-synthetic contaminated data
        % -----------------------------------------------------------------

        nTrialsUse = min(size(dataNoStimAll,1),size(dataHighStimAll,1));

        cleanTarget = dataNoStimAll(1:nTrialsUse,:);

        highStimForArtifact = dataHighStimAll(1:nTrialsUse,:);

        % Use no-stim median as a rough non-artifact reference.
        % This creates a high-stim-derived artifact-like waveform library
        % without using any artifact-removal method.
        noStimMedian = median(dataNoStimAll,1);

        artifactLibrary = highStimForArtifact - noStimMedian;

        % Randomly permute artifact trials so clean no-stim trial i is not
        % paired in a systematic way with high-stim trial i.
        permIdx = randperm(nTrialsUse);
        artifactToAdd = artifactLibrary(permIdx,:);

        semiSyntheticData = cleanTarget + artifactScale * artifactToAdd;

        %% ----------------------------------------------------------------
        % Run methods on semi-synthetic data
        % -----------------------------------------------------------------

        erpOut = ERPSubtraction(semiSyntheticData,timeVals,erpParams);
        erpCleaned = getCleanedData(erpOut);

        pcaOut = PCATemplate(semiSyntheticData,timeVals,pcaParamsK10);
        pcaCleaned = getCleanedData(pcaOut);

        smartaK3Out = SMARTALite(semiSyntheticData,timeVals,smartaK3Opts);
        smartaK3Cleaned = getCleanedData(smartaK3Out);

        smartaK5Out = SMARTALite(semiSyntheticData,timeVals,smartaK5Opts);
        smartaK5Cleaned = getCleanedData(smartaK5Out);

        ensembleCleaned = 0.5 * (smartaK3Cleaned + smartaK5Cleaned);

        methodData = { ...
            semiSyntheticData, ...
            erpCleaned, ...
            pcaCleaned, ...
            smartaK3Cleaned, ...
            smartaK5Cleaned, ...
            ensembleCleaned};

        %% ----------------------------------------------------------------
        % Reference metrics
        % -----------------------------------------------------------------

        fftCleanTarget = compute_fft_summary(cleanTarget,timeVals,fftIdx);

        freqMask = fftCleanTarget.freqAxis >= freqRangeForMetric(1) & ...
                   fftCleanTarget.freqAxis <= freqRangeForMetric(2);

        harmonicMask = buildHarmonicMask( ...
            fftCleanTarget.freqAxis, ...
            harmonicsToUse, ...
            harmonicHalfBandwidthHz, ...
            freqRangeForMetric);

        trueArtifactInWindow = semiSyntheticData(:,timeMetricIdx) - cleanTarget(:,timeMetricIdx);
        rawArtifactNorm = norm(trueArtifactInWindow,'fro');

        if rawArtifactNorm == 0
            rawArtifactNorm = eps;
        end

        cleanRMS = computeWindowRMS(cleanTarget,timeMetricIdx);
        cleanTrialSTD = computeWindowTrialSTD(cleanTarget,timeMetricIdx);

        rawFFT = compute_fft_summary(semiSyntheticData,timeVals,fftIdx);
        rawSpectralDiag = computeSpectralDiagnostics(rawFFT,fftCleanTarget,freqMask,lambdaBelow);
        rawHarmonicDiag = computeHarmonicDiagnostics(rawFFT,fftCleanTarget,harmonicMask);

        rawTotalFFTError = rawSpectralDiag.totalError;
        rawBalancedScore = rawSpectralDiag.balancedScore;
        rawHarmonicExcess = rawHarmonicDiag.harmonicAboveError;

        fprintf('%-32s normTimeErr totalFFT balanced harmonicExcess rmsErr trialSTDErr\n','Method');

        %% ----------------------------------------------------------------
        % Compute metrics for each method
        % -----------------------------------------------------------------

        for iMethod = 1:length(methodNames)

            methodName = methodNames{iMethod};
            thisData = methodData{iMethod};

            thisFFT = compute_fft_summary(thisData,timeVals,fftIdx);

            spectralDiag = computeSpectralDiagnostics(thisFFT,fftCleanTarget,freqMask,lambdaBelow);
            harmonicDiag = computeHarmonicDiagnostics(thisFFT,fftCleanTarget,harmonicMask);

            recoveryResidual = thisData(:,timeMetricIdx) - cleanTarget(:,timeMetricIdx);

            normalizedTimeRecoveryError = norm(recoveryResidual,'fro') / rawArtifactNorm;

            thisRMS = computeWindowRMS(thisData,timeMetricIdx);
            thisTrialSTD = computeWindowTrialSTD(thisData,timeMetricIdx);

            rmsErrorVsClean = abs(thisRMS - cleanRMS);
            trialSTDErrorVsClean = abs(thisTrialSTD - cleanTrialSTD);

            if rawTotalFFTError > 0
                totalFFTImprovementPct = 100 * ...
                    (rawTotalFFTError - spectralDiag.totalError) / rawTotalFFTError;
            else
                totalFFTImprovementPct = NaN;
            end

            if rawBalancedScore > 0
                balancedImprovementPct = 100 * ...
                    (rawBalancedScore - spectralDiag.balancedScore) / rawBalancedScore;
            else
                balancedImprovementPct = NaN;
            end

            if rawHarmonicExcess > 0
                harmonicSuppressionPct = 100 * ...
                    (rawHarmonicExcess - harmonicDiag.harmonicAboveError) / rawHarmonicExcess;
            else
                harmonicSuppressionPct = NaN;
            end

            rowCounter = rowCounter + 1;

            longRows(rowCounter).electrode = elecNum;
            longRows(rowCounter).method = string(methodName);

            longRows(rowCounter).normalizedTimeRecoveryError = normalizedTimeRecoveryError;

            longRows(rowCounter).totalFFTError = spectralDiag.totalError;
            longRows(rowCounter).aboveReferenceError = spectralDiag.aboveReferenceError;
            longRows(rowCounter).belowReferenceError = spectralDiag.belowReferenceError;
            longRows(rowCounter).fractionBelowReference = spectralDiag.fractionBelowReference;
            longRows(rowCounter).balancedScore = spectralDiag.balancedScore;

            longRows(rowCounter).harmonicTotalError = harmonicDiag.harmonicTotalError;
            longRows(rowCounter).harmonicAboveError = harmonicDiag.harmonicAboveError;
            longRows(rowCounter).harmonicBelowError = harmonicDiag.harmonicBelowError;

            longRows(rowCounter).rmsErrorVsClean = rmsErrorVsClean;
            longRows(rowCounter).trialSTDErrorVsClean = trialSTDErrorVsClean;

            longRows(rowCounter).totalFFTImprovementPct = totalFFTImprovementPct;
            longRows(rowCounter).balancedImprovementPct = balancedImprovementPct;
            longRows(rowCounter).harmonicSuppressionPct = harmonicSuppressionPct;

            fprintf('%-32s %.4f      %.4f   %.4f   %.4f         %.4f %.4f\n', ...
                methodName, ...
                normalizedTimeRecoveryError, ...
                spectralDiag.totalError, ...
                spectralDiag.balancedScore, ...
                harmonicDiag.harmonicAboveError, ...
                rmsErrorVsClean, ...
                trialSTDErrorVsClean);
        end

    catch ME

        warning('Failed on elec%d: %s',elecNum,ME.message);

    end
end

if isempty(longRows)
    error('No electrodes were processed successfully.');
end

longTable = struct2table(longRows);

fprintf('\nDemo 24 semi-synthetic long results table:\n');
disp(longTable);

%% ========================================================================
% Summary table
% ========================================================================

summaryRows = struct([]);
summaryCounter = 0;

for iMethod = 1:length(methodNames)

    methodName = string(methodNames{iMethod});
    idx = longTable.method == methodName;

    Tm = longTable(idx,:);

    summaryCounter = summaryCounter + 1;

    summaryRows(summaryCounter).method = methodName;
    summaryRows(summaryCounter).nElectrodes = height(Tm);

    summaryRows(summaryCounter).meanNormalizedTimeRecoveryError = mean(Tm.normalizedTimeRecoveryError);
    summaryRows(summaryCounter).stdNormalizedTimeRecoveryError = std(Tm.normalizedTimeRecoveryError);

    summaryRows(summaryCounter).meanTotalFFTError = mean(Tm.totalFFTError);
    summaryRows(summaryCounter).stdTotalFFTError = std(Tm.totalFFTError);

    summaryRows(summaryCounter).meanAboveReferenceError = mean(Tm.aboveReferenceError);
    summaryRows(summaryCounter).meanBelowReferenceError = mean(Tm.belowReferenceError);
    summaryRows(summaryCounter).meanFractionBelowReference = mean(Tm.fractionBelowReference);

    summaryRows(summaryCounter).meanBalancedScore = mean(Tm.balancedScore);
    summaryRows(summaryCounter).stdBalancedScore = std(Tm.balancedScore);

    summaryRows(summaryCounter).meanHarmonicAboveError = mean(Tm.harmonicAboveError);
    summaryRows(summaryCounter).stdHarmonicAboveError = std(Tm.harmonicAboveError);

    summaryRows(summaryCounter).meanRMSErrorVsClean = mean(Tm.rmsErrorVsClean);
    summaryRows(summaryCounter).stdRMSErrorVsClean = std(Tm.rmsErrorVsClean);

    summaryRows(summaryCounter).meanTrialSTDErrorVsClean = mean(Tm.trialSTDErrorVsClean);
    summaryRows(summaryCounter).stdTrialSTDErrorVsClean = std(Tm.trialSTDErrorVsClean);

    summaryRows(summaryCounter).meanTotalFFTImprovementPct = mean(Tm.totalFFTImprovementPct,'omitnan');
    summaryRows(summaryCounter).meanBalancedImprovementPct = mean(Tm.balancedImprovementPct,'omitnan');
    summaryRows(summaryCounter).meanHarmonicSuppressionPct = mean(Tm.harmonicSuppressionPct,'omitnan');
end

summaryTable = struct2table(summaryRows);
summaryTable = sortrows(summaryTable,'meanBalancedScore','ascend');

fprintf('\nDemo 24 semi-synthetic validation summary table:\n');
disp(summaryTable);

%% ========================================================================
% Best method counts
% ========================================================================

bestRows = struct([]);

processedElectrodes = unique(longTable.electrode);

for iElec = 1:length(processedElectrodes)

    elecNum = processedElectrodes(iElec);

    idxElec = longTable.electrode == elecNum & ismember(cellstr(longTable.method),candidateMethods);
    Te = longTable(idxElec,:);

    [~,idxTime] = min(Te.normalizedTimeRecoveryError);
    [~,idxTotal] = min(Te.totalFFTError);
    [~,idxBalanced] = min(Te.balancedScore);
    [~,idxHarmonic] = min(Te.harmonicAboveError);
    [~,idxRMS] = min(Te.rmsErrorVsClean);
    [~,idxSTD] = min(Te.trialSTDErrorVsClean);

    bestRows(iElec).electrode = elecNum;
    bestRows(iElec).bestTimeRecovery = Te.method(idxTime);
    bestRows(iElec).bestTotalFFT = Te.method(idxTotal);
    bestRows(iElec).bestBalanced = Te.method(idxBalanced);
    bestRows(iElec).bestHarmonic = Te.method(idxHarmonic);
    bestRows(iElec).bestRMS = Te.method(idxRMS);
    bestRows(iElec).bestTrialSTD = Te.method(idxSTD);
end

bestTable = struct2table(bestRows);

fprintf('\nBest method per electrode by semi-synthetic metric:\n');
disp(bestTable);

fprintf('\nBest counts by normalized time recovery error:\n');
printBestCounts(bestTable.bestTimeRecovery,candidateMethods);

fprintf('\nBest counts by total FFT recovery error:\n');
printBestCounts(bestTable.bestTotalFFT,candidateMethods);

fprintf('\nBest counts by balanced recovery score:\n');
printBestCounts(bestTable.bestBalanced,candidateMethods);

fprintf('\nBest counts by harmonic excess error:\n');
printBestCounts(bestTable.bestHarmonic,candidateMethods);

fprintf('\nBest counts by RMS error:\n');
printBestCounts(bestTable.bestRMS,candidateMethods);

fprintf('\nBest counts by trial-STD error:\n');
printBestCounts(bestTable.bestTrialSTD,candidateMethods);

%% ========================================================================
% Paired statistical comparisons
% ========================================================================

statsTable = computePairwiseStats(longTable, ...
    'SMARTALite_Ensemble_K3K5', ...
    {'ERPSubtraction','PCATemplate_K10','SMARTALite_K3','SMARTALite_K5'}, ...
    {'normalizedTimeRecoveryError','totalFFTError','balancedScore','harmonicAboveError','rmsErrorVsClean','trialSTDErrorVsClean'});

fprintf('\nPaired statistics: Ensemble K3K5 vs comparison methods\n');
disp(statsTable);

%% ========================================================================
% Figures
% ========================================================================

candidateSummary = summaryTable(~strcmp(cellstr(summaryTable.method),'RawSemiSynthetic'),:);

figure('Name','demo24 semi-synthetic recovery errors');

barData = [ ...
    candidateSummary.meanNormalizedTimeRecoveryError, ...
    candidateSummary.meanTotalFFTError, ...
    candidateSummary.meanBalancedScore];

bar(categorical(cellstr(candidateSummary.method)),barData);
ylabel('Metric value');
title('Semi-synthetic validation: recovery errors');
legend({'Normalized time recovery error','Total FFT error','Balanced score'}, ...
    'Location','best');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo24_recovery_errors.png'));
    savefig(gcf,fullfile(figFolder,'demo24_recovery_errors.fig'));
end

figure('Name','demo24 semi-synthetic over-cleaning');

barData = [ ...
    candidateSummary.meanAboveReferenceError, ...
    candidateSummary.meanBelowReferenceError];

bar(categorical(cellstr(candidateSummary.method)),barData);
ylabel('Spectral error norm');
title('Semi-synthetic validation: above vs below clean reference');
legend({'Above clean reference','Below clean reference'}, ...
    'Location','best');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo24_above_below_clean_reference.png'));
    savefig(gcf,fullfile(figFolder,'demo24_above_below_clean_reference.fig'));
end

figure('Name','demo24 harmonic suppression');

bar(categorical(cellstr(candidateSummary.method)),candidateSummary.meanHarmonicSuppressionPct);
ylabel('Harmonic suppression relative to raw semi-synthetic (%)');
title('Semi-synthetic validation: stimulation harmonic suppression');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo24_harmonic_suppression.png'));
    savefig(gcf,fullfile(figFolder,'demo24_harmonic_suppression.fig'));
end

figure('Name','demo24 time-domain recovery');

barData = [ ...
    candidateSummary.meanRMSErrorVsClean, ...
    candidateSummary.meanTrialSTDErrorVsClean];

bar(categorical(cellstr(candidateSummary.method)),barData);
ylabel('Metric value');
title('Semi-synthetic validation: time-domain recovery');
legend({'RMS error vs clean','Trial-STD error vs clean'}, ...
    'Location','best');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo24_time_domain_recovery.png'));
    savefig(gcf,fullfile(figFolder,'demo24_time_domain_recovery.fig'));
end

%% ========================================================================
% Save results
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo24_semisynthetic_validation_' timestampString '.mat']);
longCsvFile = fullfile(resultsFolder,['demo24_semisynthetic_validation_long_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo24_semisynthetic_validation_summary_' timestampString '.csv']);
bestCsvFile = fullfile(resultsFolder,['demo24_semisynthetic_validation_best_' timestampString '.csv']);
statsCsvFile = fullfile(resultsFolder,['demo24_semisynthetic_validation_stats_' timestampString '.csv']);

save(resultsMatFile, ...
    'longTable', ...
    'summaryTable', ...
    'bestTable', ...
    'statsTable', ...
    'methodNames', ...
    'candidateMethods', ...
    'goodV1Electrodes', ...
    'pulseTimes', ...
    'smartaStimFreq', ...
    'smartaPrePulse', ...
    'artifactWindow', ...
    'fftWindow', ...
    'timeMetricWindow', ...
    'freqRangeForMetric', ...
    'harmonicsToUse', ...
    'harmonicHalfBandwidthHz', ...
    'artifactScale', ...
    'lambdaBelow', ...
    'rngSeed');

writetable(longTable,longCsvFile);
writetable(summaryTable,summaryCsvFile);
writetable(bestTable,bestCsvFile);
writetable(statsTable,statsCsvFile);

fprintf('\nSaved demo_24 results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',longCsvFile);
fprintf('%s\n',summaryCsvFile);
fprintf('%s\n',bestCsvFile);
fprintf('%s\n',statsCsvFile);

if saveFigures
    fprintf('\nSaved demo_24 figures in:\n');
    fprintf('%s\n',figFolder);
end

fprintf('\n============================================================\n');
fprintf('demo_24 semi-synthetic validation complete\n');
fprintf('============================================================\n');

fprintf('\nInterpretation guide:\n');
fprintf('If SMARTALite Ensemble K3K5 has the lowest total FFT or balanced score,\n');
fprintf('then the final method is supported even when the clean target is known.\n');
fprintf('If SMARTALite K3 wins the harmonic metric, that is consistent with demo_22.\n');

%% ========================================================================
% Local helper functions
% ========================================================================

function cleanedData = getCleanedData(methodOut)

    if isfield(methodOut,'cleanedData')
        cleanedData = methodOut.cleanedData;
    elseif isfield(methodOut,'cleanedTrials')
        cleanedData = methodOut.cleanedTrials;
    else
        error('Method output has neither cleanedData nor cleanedTrials.');
    end
end

function diagOut = computeSpectralDiagnostics(fftCleaned,fftCleanTarget,freqMask,lambdaBelow)

    diffVals = fftCleaned.logMeanMagnitude(freqMask) - ...
               fftCleanTarget.logMeanMagnitude(freqMask);

    aboveVals = max(diffVals,0);
    belowVals = max(-diffVals,0);

    diagOut = struct();

    diagOut.totalError = norm(diffVals);
    diagOut.aboveReferenceError = norm(aboveVals);
    diagOut.belowReferenceError = norm(belowVals);
    diagOut.fractionBelowReference = mean(diffVals < 0);
    diagOut.meanSignedDiff = mean(diffVals);

    diagOut.balancedScore = diagOut.aboveReferenceError + ...
                            lambdaBelow * diagOut.belowReferenceError;
end

function harmonicMask = buildHarmonicMask(freqAxis,harmonicsToUse,halfBandwidthHz,freqRangeForMetric)

    harmonicMask = false(size(freqAxis));

    for iH = 1:length(harmonicsToUse)

        h = harmonicsToUse(iH);

        if h < freqRangeForMetric(1) || h > freqRangeForMetric(2)
            continue;
        end

        harmonicMask = harmonicMask | abs(freqAxis - h) <= halfBandwidthHz;
    end
end

function diagOut = computeHarmonicDiagnostics(fftCleaned,fftCleanTarget,harmonicMask)

    diffVals = fftCleaned.logMeanMagnitude(harmonicMask) - ...
               fftCleanTarget.logMeanMagnitude(harmonicMask);

    aboveVals = max(diffVals,0);
    belowVals = max(-diffVals,0);

    diagOut = struct();

    diagOut.harmonicTotalError = norm(diffVals);
    diagOut.harmonicAboveError = norm(aboveVals);
    diagOut.harmonicBelowError = norm(belowVals);
    diagOut.harmonicFractionBelow = mean(diffVals < 0);
end

function rmsVal = computeWindowRMS(data,idx)

    x = data(:,idx);
    rmsVal = sqrt(mean(x(:).^2));
end

function trialSTDVal = computeWindowTrialSTD(data,idx)

    x = data(:,idx);
    stdTime = std(x,0,1);
    trialSTDVal = mean(stdTime);
end

function printBestCounts(bestMethods,candidateMethods)

    bestMethods = string(bestMethods);

    for i = 1:length(candidateMethods)
        methodName = string(candidateMethods{i});
        nBest = sum(bestMethods == methodName);
        fprintf('%s: %d electrodes\n',candidateMethods{i},nBest);
    end
end

function statsTable = computePairwiseStats(longTable,mainMethod,comparisonMethods,metricNames)

    rows = struct([]);
    r = 0;

    for iComp = 1:length(comparisonMethods)

        compMethod = comparisonMethods{iComp};

        for iMetric = 1:length(metricNames)

            metricName = metricNames{iMetric};

            mainVals = getMetricByMethod(longTable,mainMethod,metricName);
            compVals = getMetricByMethod(longTable,compMethod,metricName);

            validIdx = ~isnan(mainVals) & ~isnan(compVals);
            mainVals = mainVals(validIdx);
            compVals = compVals(validIdx);

            diffVals = compVals - mainVals;
            % Positive diff means main method has smaller metric value.

            r = r + 1;

            rows(r).mainMethod = string(mainMethod);
            rows(r).comparisonMethod = string(compMethod);
            rows(r).metric = string(metricName);

            rows(r).n = length(diffVals);
            rows(r).meanComparisonMinusMain = mean(diffVals);
            rows(r).medianComparisonMinusMain = median(diffVals);
            rows(r).nMainBetter = sum(diffVals > 0);
            rows(r).nComparisonBetter = sum(diffVals < 0);
            rows(r).percentMainBetter = 100 * sum(diffVals > 0) / length(diffVals);

            try
                pVal = signrank(mainVals,compVals);
            catch
                pVal = NaN;
            end

            rows(r).signrankP = pVal;
        end
    end

    statsTable = struct2table(rows);
end

function vals = getMetricByMethod(longTable,methodName,metricName)

    idx = longTable.method == string(methodName);
    T = longTable(idx,:);

    [~,sortIdx] = sort(T.electrode);
    T = T(sortIdx,:);

    vals = T.(metricName);
end