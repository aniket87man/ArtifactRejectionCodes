clear; close all; clc;

%% demo_22_FinalValidationMetrics
%
% Goal:
% Final quantitative validation of ICMS artifact-removal methods.
%
% Methods:
%   1. Raw high-stim
%   2. ERPSubtraction
%   3. PCATemplate K=20
%   4. SMARTALite K=3
%   5. SMARTALite K=5
%   6. SMARTALite Ensemble K3K5
%
% Metrics:
%   1. Total FFT error vs no-stim reference
%   2. Above-reference error: residual excess spectral power
%   3. Below-reference error: possible over-cleaning
%   4. Balanced score: above + lambda * below
%   5. Harmonic excess error at 20 Hz stimulation harmonics
%   6. Harmonic suppression relative to raw high-stim
%   7. Time-domain RMS error vs no-stim
%   8. Trial-variability error vs no-stim
%   9. No-stim distortion control
%
% Important:
% no-stim is a reference condition, not absolute ground truth.

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

% FFT metric settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% Stimulation harmonic settings
stimFreq = 20;
harmonicsToUse = stimFreq:stimFreq:200;
harmonicHalfBandwidthHz = 1.0;

% Time-domain metric window
timeMetricWindow = [0 0.4];

% PCA settings
artifactWindow = [0 0.4];

% SMARTALite full-cycle settings
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Over-cleaning penalty
lambdaBelow = 1.5;

saveFigures = true;

%% ========================================================================
% Paths and loading
% ========================================================================

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

fftIdx = find(timeVals > fftWindow(1) & timeVals <= fftWindow(2));
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
% Corrected good V1 electrode selection
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

fprintf('\nFinal good V1 electrodes used in demo_22:\n');
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
erpParams.subtractWindow = artifactWindow;
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
    'RawHighStim', ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_Ensemble_K3K5'};

%% ========================================================================
% Output folders
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo22_FinalValidationMetrics');

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

        dataNoStim = analogData(noStimTrials,:);
        dataHighStim = analogData(highStimTrials,:);

        %% Clean high-stim data

        erpOutHigh = ERPSubtraction(dataHighStim,timeVals,erpParams);
        erpHigh = getCleanedData(erpOutHigh);

        pcaOutHigh = PCATemplate(dataHighStim,timeVals,pcaParamsK10);
        pcaHigh = getCleanedData(pcaOutHigh);

        smartaK3OutHigh = SMARTALite(dataHighStim,timeVals,smartaK3Opts);
        smartaK3High = getCleanedData(smartaK3OutHigh);

        smartaK5OutHigh = SMARTALite(dataHighStim,timeVals,smartaK5Opts);
        smartaK5High = getCleanedData(smartaK5OutHigh);

        ensembleHigh = 0.5 * (smartaK3High + smartaK5High);

        methodData = { ...
            dataHighStim, ...
            erpHigh, ...
            pcaHigh, ...
            smartaK3High, ...
            smartaK5High, ...
            ensembleHigh};

        %% No-stim distortion control
        % Apply each method to no-stim trials and measure how much no-stim changes.
        % This is a sensitivity/control metric, not the primary selection metric.

        noStimOriginal = dataNoStim;

        erpOutNoStim = ERPSubtraction(dataNoStim,timeVals,erpParams);
        erpNoStimCleaned = getCleanedData(erpOutNoStim);

        pcaOutNoStim = PCATemplate(dataNoStim,timeVals,pcaParamsK10);
        pcaNoStimCleaned = getCleanedData(pcaOutNoStim);

        smartaK3OutNoStim = SMARTALite(dataNoStim,timeVals,smartaK3Opts);
        smartaK3NoStimCleaned = getCleanedData(smartaK3OutNoStim);

        smartaK5OutNoStim = SMARTALite(dataNoStim,timeVals,smartaK5Opts);
        smartaK5NoStimCleaned = getCleanedData(smartaK5OutNoStim);

        ensembleNoStimCleaned = 0.5 * (smartaK3NoStimCleaned + smartaK5NoStimCleaned);

        noStimDistortionData = { ...
            noStimOriginal, ...
            erpNoStimCleaned, ...
            pcaNoStimCleaned, ...
            smartaK3NoStimCleaned, ...
            smartaK5NoStimCleaned, ...
            ensembleNoStimCleaned};

        %% Reference FFT and baseline metrics

        fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);

        freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
                   fftNoStim.freqAxis <= freqRangeForMetric(2);

        harmonicMask = buildHarmonicMask( ...
            fftNoStim.freqAxis, ...
            harmonicsToUse, ...
            harmonicHalfBandwidthHz, ...
            freqRangeForMetric);

        noStimRMS = computeWindowRMS(dataNoStim,timeMetricIdx);
        noStimTrialSTD = computeWindowTrialSTD(dataNoStim,timeMetricIdx);

        rawFFT = compute_fft_summary(dataHighStim,timeVals,fftIdx);
        rawHarmonicDiag = computeHarmonicDiagnostics(rawFFT,fftNoStim,harmonicMask);
        rawHarmonicExcess = rawHarmonicDiag.harmonicAboveError;

        fprintf('%-30s totalErr balanced harmonicExcess harmonicSupp RMSerr STDerr noStimDist\n','Method');

        %% Compute metrics for each method

        for iMethod = 1:length(methodNames)

            methodName = methodNames{iMethod};
            thisData = methodData{iMethod};

            thisFFT = compute_fft_summary(thisData,timeVals,fftIdx);

            spectralDiag = computeSpectralDiagnostics(thisFFT,fftNoStim,freqMask,lambdaBelow);
            harmonicDiag = computeHarmonicDiagnostics(thisFFT,fftNoStim,harmonicMask);

            if rawHarmonicExcess > 0
                harmonicSuppressionPct = 100 * ...
                    (rawHarmonicExcess - harmonicDiag.harmonicAboveError) / rawHarmonicExcess;
            else
                harmonicSuppressionPct = NaN;
            end

            thisRMS = computeWindowRMS(thisData,timeMetricIdx);
            thisTrialSTD = computeWindowTrialSTD(thisData,timeMetricIdx);

            rmsError = abs(thisRMS - noStimRMS);
            trialSTDError = abs(thisTrialSTD - noStimTrialSTD);

            noStimDistortion = computeNoStimDistortion( ...
                noStimOriginal, ...
                noStimDistortionData{iMethod}, ...
                timeVals, ...
                fftIdx, ...
                freqMask);

            rowCounter = rowCounter + 1;

            longRows(rowCounter).electrode = elecNum;
            longRows(rowCounter).method = string(methodName);

            longRows(rowCounter).totalFFTError = spectralDiag.totalError;
            longRows(rowCounter).aboveReferenceError = spectralDiag.aboveReferenceError;
            longRows(rowCounter).belowReferenceError = spectralDiag.belowReferenceError;
            longRows(rowCounter).fractionBelowReference = spectralDiag.fractionBelowReference;
            longRows(rowCounter).balancedScore = spectralDiag.balancedScore;

            longRows(rowCounter).harmonicTotalError = harmonicDiag.harmonicTotalError;
            longRows(rowCounter).harmonicAboveError = harmonicDiag.harmonicAboveError;
            longRows(rowCounter).harmonicBelowError = harmonicDiag.harmonicBelowError;
            longRows(rowCounter).harmonicSuppressionPct = harmonicSuppressionPct;

            longRows(rowCounter).windowRMS = thisRMS;
            longRows(rowCounter).rmsErrorVsNoStim = rmsError;

            longRows(rowCounter).trialSTD = thisTrialSTD;
            longRows(rowCounter).trialSTDErrorVsNoStim = trialSTDError;

            longRows(rowCounter).noStimDistortion = noStimDistortion;

            fprintf('%-30s %.4f   %.4f   %.4f        %.2f       %.4f %.4f %.4f\n', ...
                methodName, ...
                spectralDiag.totalError, ...
                spectralDiag.balancedScore, ...
                harmonicDiag.harmonicAboveError, ...
                harmonicSuppressionPct, ...
                rmsError, ...
                trialSTDError, ...
                noStimDistortion);
        end

    catch ME

        warning('Failed on elec%d: %s',elecNum,ME.message);

    end
end

if isempty(longRows)
    error('No electrodes were processed successfully.');
end

longTable = struct2table(longRows);

fprintf('\nDemo 22 long results table:\n');
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

    summaryRows(summaryCounter).meanTotalFFTError = mean(Tm.totalFFTError);
    summaryRows(summaryCounter).stdTotalFFTError = std(Tm.totalFFTError);

    summaryRows(summaryCounter).meanAboveReferenceError = mean(Tm.aboveReferenceError);
    summaryRows(summaryCounter).meanBelowReferenceError = mean(Tm.belowReferenceError);
    summaryRows(summaryCounter).meanFractionBelowReference = mean(Tm.fractionBelowReference);

    summaryRows(summaryCounter).meanBalancedScore = mean(Tm.balancedScore);
    summaryRows(summaryCounter).stdBalancedScore = std(Tm.balancedScore);

    summaryRows(summaryCounter).meanHarmonicAboveError = mean(Tm.harmonicAboveError);
    summaryRows(summaryCounter).stdHarmonicAboveError = std(Tm.harmonicAboveError);

    summaryRows(summaryCounter).meanHarmonicSuppressionPct = mean(Tm.harmonicSuppressionPct,'omitnan');
    summaryRows(summaryCounter).stdHarmonicSuppressionPct = std(Tm.harmonicSuppressionPct,'omitnan');

    summaryRows(summaryCounter).meanRMSErrorVsNoStim = mean(Tm.rmsErrorVsNoStim);
    summaryRows(summaryCounter).stdRMSErrorVsNoStim = std(Tm.rmsErrorVsNoStim);

    summaryRows(summaryCounter).meanTrialSTDErrorVsNoStim = mean(Tm.trialSTDErrorVsNoStim);
    summaryRows(summaryCounter).stdTrialSTDErrorVsNoStim = std(Tm.trialSTDErrorVsNoStim);

    summaryRows(summaryCounter).meanNoStimDistortion = mean(Tm.noStimDistortion);
    summaryRows(summaryCounter).stdNoStimDistortion = std(Tm.noStimDistortion);
end

summaryTable = struct2table(summaryRows);
summaryTable = sortrows(summaryTable,'meanBalancedScore','ascend');

fprintf('\nDemo 22 final validation summary table:\n');
disp(summaryTable);

%% ========================================================================
% Best method counts
% ========================================================================

candidateMethods = { ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_Ensemble_K3K5'};

bestRows = struct([]);

for iElec = 1:length(elecNums)

    elecNum = elecNums(iElec);

    idxElec = longTable.electrode == elecNum & ismember(cellstr(longTable.method),candidateMethods);
    Te = longTable(idxElec,:);

    [~,idxTotal] = min(Te.totalFFTError);
    [~,idxBalanced] = min(Te.balancedScore);
    [~,idxHarmonic] = min(Te.harmonicAboveError);
    [~,idxRMS] = min(Te.rmsErrorVsNoStim);
    [~,idxSTD] = min(Te.trialSTDErrorVsNoStim);

    bestRows(iElec).electrode = elecNum;
    bestRows(iElec).bestTotalFFT = Te.method(idxTotal);
    bestRows(iElec).bestBalanced = Te.method(idxBalanced);
    bestRows(iElec).bestHarmonicSuppression = Te.method(idxHarmonic);
    bestRows(iElec).bestRMSError = Te.method(idxRMS);
    bestRows(iElec).bestTrialSTDError = Te.method(idxSTD);
end

bestTable = struct2table(bestRows);

fprintf('\nBest method per electrode by metric:\n');
disp(bestTable);

fprintf('\nBest counts by total FFT error:\n');
printBestCounts(bestTable.bestTotalFFT,candidateMethods);

fprintf('\nBest counts by balanced over-cleaning score:\n');
printBestCounts(bestTable.bestBalanced,candidateMethods);

fprintf('\nBest counts by harmonic excess error:\n');
printBestCounts(bestTable.bestHarmonicSuppression,candidateMethods);

fprintf('\nBest counts by RMS error vs no-stim:\n');
printBestCounts(bestTable.bestRMSError,candidateMethods);

fprintf('\nBest counts by trial-STD error vs no-stim:\n');
printBestCounts(bestTable.bestTrialSTDError,candidateMethods);

%% ========================================================================
% Paired statistical comparisons
% ========================================================================

statsTable = computePairwiseStats(longTable, ...
    'SMARTALite_Ensemble_K3K5', ...
    {'PCATemplate_K10','SMARTALite_K3','SMARTALite_K5'}, ...
    {'totalFFTError','balancedScore','harmonicAboveError','rmsErrorVsNoStim','trialSTDErrorVsNoStim'});

fprintf('\nPaired statistics: Ensemble K3K5 vs comparison methods\n');
disp(statsTable);

%% ========================================================================
% Figures
% ========================================================================

candidateSummary = summaryTable(~strcmp(cellstr(summaryTable.method),'RawHighStim'),:);

figure('Name','demo22 final validation: spectral metrics');

barData = [ ...
    candidateSummary.meanTotalFFTError, ...
    candidateSummary.meanBalancedScore, ...
    candidateSummary.meanHarmonicAboveError];

bar(categorical(cellstr(candidateSummary.method)),barData);
ylabel('Metric value');
title('Final validation: spectral metrics');
legend({'Total FFT error','Balanced score','Harmonic excess error'},'Location','best');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo22_spectral_metrics.png'));
    savefig(gcf,fullfile(figFolder,'demo22_spectral_metrics.fig'));
end

figure('Name','demo22 final validation: over-cleaning');

barData = [ ...
    candidateSummary.meanAboveReferenceError, ...
    candidateSummary.meanBelowReferenceError];

bar(categorical(cellstr(candidateSummary.method)),barData);
ylabel('Norm of spectral difference');
title('Final validation: above-reference vs below-reference error');
legend({'Above reference / residual excess','Below reference / possible over-cleaning'}, ...
    'Location','best');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo22_overcleaning_metrics.png'));
    savefig(gcf,fullfile(figFolder,'demo22_overcleaning_metrics.fig'));
end

figure('Name','demo22 final validation: time-domain metrics');

barData = [ ...
    candidateSummary.meanRMSErrorVsNoStim, ...
    candidateSummary.meanTrialSTDErrorVsNoStim, ...
    candidateSummary.meanNoStimDistortion];

bar(categorical(cellstr(candidateSummary.method)),barData);
ylabel('Metric value');
title('Final validation: time-domain and no-stim control metrics');
legend({'RMS error vs no-stim','Trial-STD error vs no-stim','No-stim distortion'}, ...
    'Location','best');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo22_time_domain_metrics.png'));
    savefig(gcf,fullfile(figFolder,'demo22_time_domain_metrics.fig'));
end

%% ========================================================================
% Save results
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo22_final_validation_metrics_' timestampString '.mat']);
longCsvFile = fullfile(resultsFolder,['demo22_final_validation_long_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo22_final_validation_summary_' timestampString '.csv']);
bestCsvFile = fullfile(resultsFolder,['demo22_final_validation_best_counts_' timestampString '.csv']);
statsCsvFile = fullfile(resultsFolder,['demo22_final_validation_stats_' timestampString '.csv']);

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
    'fftWindow', ...
    'freqRangeForMetric', ...
    'timeMetricWindow', ...
    'harmonicsToUse', ...
    'harmonicHalfBandwidthHz', ...
    'lambdaBelow');

writetable(longTable,longCsvFile);
writetable(summaryTable,summaryCsvFile);
writetable(bestTable,bestCsvFile);
writetable(statsTable,statsCsvFile);

fprintf('\nSaved demo_22 results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',longCsvFile);
fprintf('%s\n',summaryCsvFile);
fprintf('%s\n',bestCsvFile);
fprintf('%s\n',statsCsvFile);

if saveFigures
    fprintf('\nSaved demo_22 figures in:\n');
    fprintf('%s\n',figFolder);
end

fprintf('\n============================================================\n');
fprintf('demo_22 final validation metrics complete\n');
fprintf('============================================================\n');

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

function diagOut = computeSpectralDiagnostics(fftCleaned,fftNoStim,freqMask,lambdaBelow)

    diffVals = fftCleaned.logMeanMagnitude(freqMask) - ...
               fftNoStim.logMeanMagnitude(freqMask);

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

function diagOut = computeHarmonicDiagnostics(fftCleaned,fftNoStim,harmonicMask)

    diffVals = fftCleaned.logMeanMagnitude(harmonicMask) - ...
               fftNoStim.logMeanMagnitude(harmonicMask);

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

function distortionVal = computeNoStimDistortion(originalNoStim,cleanedNoStim,timeVals,fftIdx,freqMask)

    fftOriginal = compute_fft_summary(originalNoStim,timeVals,fftIdx);
    fftCleaned = compute_fft_summary(cleanedNoStim,timeVals,fftIdx);

    diffVals = fftCleaned.logMeanMagnitude(freqMask) - ...
               fftOriginal.logMeanMagnitude(freqMask);

    distortionVal = norm(diffVals);
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