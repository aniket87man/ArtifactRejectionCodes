clear; close all; clc;

%% demo_16_SMARTA
%
% Goal:
% Test SMARTA-style methods on ICMS artifact removal.
%
% Methods compared:
% 1. Raw high-stim
% 2. ERPSubtraction
% 3. ERPAligned
% 4. TunedHybrid
% 5. PCATemplate K=3
% 6. PCATemplate K=20
% 7. SMARTALite
% 8. SMARTAFull
%
% SMARTA paper idea:
% segment artifacts, find similar artifact segments using a denoised
% similarity metric, build a target-specific median template, subtract it.
%
% Our adaptation:
% pulse-wise SMARTA for ICMS trial data.

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

% First SMARTA test electrodes
% elec1 = stimulation electrode
% elec3,6,35 = strong PCA K=20 examples
% elec8,24 = exceptions from demo_15
electrodesToTest = [1 3 6 8 24 35];

% Artifact/FFT settings
artifactWindow = [-0.02 0.4];
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% Approximate pulse times for high-stim condition
% Adjust later if needed after visual inspection.

pulseTimes = 0 + (0:7)*0.05;   % 8 pulse/cycle starts: 0, 0.05, ..., 0.35 s

smartaK = 5;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Save figures?
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

fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));
fprintf('FFT window samples = %d\n',length(fftIdx));
fprintf('FFT window duration = %.4f s\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% ========================================================================
% Confirm required methods are visible
% ========================================================================

requiredMethods = { ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'ERPAlignedPulsewise', ...
    'PCATemplate', ...
    'SMARTALite', ...
    'SMARTAFull', ...
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
% Method parameters
% ========================================================================

% ERPSubtraction
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

% ERPAligned
alignedParams = struct();
alignedParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
alignedParams.alignWindow = [-0.01 0.03];
alignedParams.maxShiftMS = 10;
alignedParams.doBaselineCorrection = false;

% Tuned hybrid from demo_14/demo_15
hybridPulseParams = struct();
hybridPulseParams.pulseSearchWindow = [-0.02 0.32];
hybridPulseParams.pulseWindowMS = [-10 40];
hybridPulseParams.minPulseDistanceMS = 25;
hybridPulseParams.expectedNumPulses = 7;
hybridPulseParams.maxNumPulses = 20;
hybridPulseParams.thresholdMAD = 2;
hybridPulseParams.localBaselineMS = 20;
hybridPulseParams.templateStatistic = 'median';
hybridPulseParams.taperEdgeMS = 2;
hybridPulseParams.doBaselineCorrection = false;

hybridParams = struct();
hybridParams.alignedParams = alignedParams;
hybridParams.pulseParams = hybridPulseParams;

% PCA K=3
pcaParamsK3 = struct();
pcaParamsK3.artifactWindow = artifactWindow;
pcaParamsK3.numComponents = 3;
pcaParamsK3.removeMeanTemplate = true;
pcaParamsK3.taperEdgeMS = 2;
pcaParamsK3.doBaselineCorrection = false;

% PCA K=20
pcaParamsK20 = struct();
pcaParamsK20.artifactWindow = artifactWindow;
pcaParamsK20.numComponents = 20;
pcaParamsK20.removeMeanTemplate = true;
pcaParamsK20.taperEdgeMS = 2;
pcaParamsK20.doBaselineCorrection = false;

% SMARTALite
smartaLiteOpts = struct();
smartaLiteOpts.pulseTimes = pulseTimes;
smartaLiteOpts.stiFreq = smartaStimFreq;
smartaLiteOpts.prePulse = smartaPrePulse;
smartaLiteOpts.K = smartaK;
smartaLiteOpts.window = fftWindow;
smartaLiteOpts.computeFFT = false;

% SMARTAFull
smartaFullOpts = smartaLiteOpts;
smartaFullOpts.hpCutoff = 150;
smartaFullOpts.shrinkLoss = 'fro';
smartaFullOpts.shrinkKL = 10;
smartaFullOpts.shrinkKH = 15;
%% ========================================================================
% Output folders
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo16_SMARTA');

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Loop over electrodes
% ========================================================================

methodNames = { ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'TunedHybrid', ...
    'PCATemplate_K3', ...
    'PCATemplate_K20', ...
    'SMARTALite', ...
    'SMARTAFull'};

resultRows = struct([]);
rowCounter = 0;

for iElec = 1:length(electrodesToTest)

    elecNum = electrodesToTest(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    fprintf('\n============================================================\n');
    fprintf('Processing elec%d\n',elecNum);
    fprintf('============================================================\n');

    if ~exist(elecFile,'file')
        warning('Could not find %s. Skipping elec%d.',elecFile,elecNum);
        continue;
    end

    D = load(elecFile);

    if ~isfield(D,'analogData')
        warning('analogData not found in elec%d. Skipping.',elecNum);
        continue;
    end

    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    %% Reference FFT
    fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
    fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);

    freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
               fftNoStim.freqAxis <= freqRangeForMetric(2);

    rawError = norm(fftHighRaw.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

    %% Existing methods

    erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);
    fftERP = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);
    erpError = norm(fftERP.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

    alignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);
    fftAligned = compute_fft_summary(alignedOut.cleanedData,timeVals,fftIdx);
    alignedError = norm(fftAligned.logMeanMagnitude(freqMask) - ...
                        fftNoStim.logMeanMagnitude(freqMask));

    hybridOut = ERPAlignedPulsewise(dataHighStim,timeVals,hybridParams);
    fftHybrid = compute_fft_summary(hybridOut.cleanedData,timeVals,fftIdx);
    hybridError = norm(fftHybrid.logMeanMagnitude(freqMask) - ...
                       fftNoStim.logMeanMagnitude(freqMask));

    pcaK3Out = PCATemplate(dataHighStim,timeVals,pcaParamsK3);
    fftPCAK3 = compute_fft_summary(pcaK3Out.cleanedData,timeVals,fftIdx);
    pcaK3Error = norm(fftPCAK3.logMeanMagnitude(freqMask) - ...
                      fftNoStim.logMeanMagnitude(freqMask));

    pcaK20Out = PCATemplate(dataHighStim,timeVals,pcaParamsK20);
    fftPCAK20 = compute_fft_summary(pcaK20Out.cleanedData,timeVals,fftIdx);
    pcaK20Error = norm(fftPCAK20.logMeanMagnitude(freqMask) - ...
                       fftNoStim.logMeanMagnitude(freqMask));

    %% SMARTALite

    smartaLiteError = NaN;
    smartaLiteImprovement = NaN;
    smartaLiteOut = [];

    try
        smartaLiteOut = SMARTALite(dataHighStim,timeVals,smartaLiteOpts);
        fftSMARTALite = compute_fft_summary(smartaLiteOut.cleanedData,timeVals,fftIdx);

        smartaLiteError = norm(fftSMARTALite.logMeanMagnitude(freqMask) - ...
                               fftNoStim.logMeanMagnitude(freqMask));
    catch ME
        warning('SMARTALite failed on elec%d: %s',elecNum,ME.message);
        smartaLiteOut = [];
        clear fftSMARTALite
    end

    %% SMARTAFull

    smartaFullError = NaN;
    smartaFullImprovement = NaN;
    smartaFullOut = [];

    try
        smartaFullOut = SMARTAFull(dataHighStim,timeVals,smartaFullOpts);
        fftSMARTAFull = compute_fft_summary(smartaFullOut.cleanedData,timeVals,fftIdx);

        smartaFullError = norm(fftSMARTAFull.logMeanMagnitude(freqMask) - ...
                               fftNoStim.logMeanMagnitude(freqMask));
    catch ME
        warning('SMARTAFull failed on elec%d: %s',elecNum,ME.message);
        smartaFullOut = [];
        clear fftSMARTAFull
    end

    %% Improvements

    erpImprovement = 100 * (rawError - erpError) / rawError;
    alignedImprovement = 100 * (rawError - alignedError) / rawError;
    hybridImprovement = 100 * (rawError - hybridError) / rawError;
    pcaK3Improvement = 100 * (rawError - pcaK3Error) / rawError;
    pcaK20Improvement = 100 * (rawError - pcaK20Error) / rawError;

    if ~isnan(smartaLiteError)
        smartaLiteImprovement = 100 * (rawError - smartaLiteError) / rawError;
    end

    if ~isnan(smartaFullError)
        smartaFullImprovement = 100 * (rawError - smartaFullError) / rawError;
    end

    %% Store row

    rowCounter = rowCounter + 1;

    resultRows(rowCounter).electrode = elecNum;
    resultRows(rowCounter).rawError = rawError;

    resultRows(rowCounter).erpError = erpError;
    resultRows(rowCounter).erpImprovement = erpImprovement;

    resultRows(rowCounter).alignedError = alignedError;
    resultRows(rowCounter).alignedImprovement = alignedImprovement;

    resultRows(rowCounter).hybridError = hybridError;
    resultRows(rowCounter).hybridImprovement = hybridImprovement;

    resultRows(rowCounter).pcaK3Error = pcaK3Error;
    resultRows(rowCounter).pcaK3Improvement = pcaK3Improvement;

    resultRows(rowCounter).pcaK20Error = pcaK20Error;
    resultRows(rowCounter).pcaK20Improvement = pcaK20Improvement;

    resultRows(rowCounter).smartaLiteError = smartaLiteError;
    resultRows(rowCounter).smartaLiteImprovement = smartaLiteImprovement;

    resultRows(rowCounter).smartaFullError = smartaFullError;
    resultRows(rowCounter).smartaFullImprovement = smartaFullImprovement;

    %% Print

    fprintf('Raw error        = %.4f\n',rawError);
    fprintf('ERPSubtraction   = %.4f, improvement %.2f %%\n',erpError,erpImprovement);
    fprintf('ERPAligned       = %.4f, improvement %.2f %%\n',alignedError,alignedImprovement);
    fprintf('TunedHybrid      = %.4f, improvement %.2f %%\n',hybridError,hybridImprovement);
    fprintf('PCATemplate K=3  = %.4f, improvement %.2f %%\n',pcaK3Error,pcaK3Improvement);
    fprintf('PCATemplate K=20 = %.4f, improvement %.2f %%\n',pcaK20Error,pcaK20Improvement);
    fprintf('SMARTALite       = %.4f, improvement %.2f %%\n',smartaLiteError,smartaLiteImprovement);
    fprintf('SMARTAFull       = %.4f, improvement %.2f %%\n',smartaFullError,smartaFullImprovement);

    %% Figure: FFT comparison

    if saveFigures

        figFFT = figure('Name',['demo16 SMARTA elec' num2str(elecNum) ' FFT']);

        plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.4);
        hold on;
        plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'Color',[0.7 0.7 0.7],'LineWidth',1.0);
        plot(fftERP.freqAxis,fftERP.logMeanMagnitude,'LineWidth',1.0);
        plot(fftPCAK20.freqAxis,fftPCAK20.logMeanMagnitude,'LineWidth',1.2);

        legendEntries = {'No-stim','Raw high-stim','ERPSubtraction','PCATemplate K=20'};

        if exist('fftSMARTALite','var') && ~isempty(smartaLiteOut)
            plot(fftSMARTALite.freqAxis,fftSMARTALite.logMeanMagnitude,'LineWidth',1.2);
            legendEntries{end+1} = 'SMARTALite';
        end

        if exist('fftSMARTAFull','var') && ~isempty(smartaFullOut)
            plot(fftSMARTAFull.freqAxis,fftSMARTAFull.logMeanMagnitude,'LineWidth',1.2);
            legendEntries{end+1} = 'SMARTAFull';
        end

        xlim(freqRangeForMetric);
        xlabel('Frequency (Hz)');
        ylabel('log10 mean FFT magnitude');
        title(['demo16 SMARTA FFT comparison, elec' num2str(elecNum)]);
        legend(legendEntries,'Location','best');
        grid on;

        saveas(figFFT,fullfile(figFolder,['demo16_elec' num2str(elecNum) '_fft.png']));
        savefig(figFFT,fullfile(figFolder,['demo16_elec' num2str(elecNum) '_fft.fig']));

        %% Figure: time series comparison

        figTS = figure('Name',['demo16 SMARTA elec' num2str(elecNum) ' time series']);

        plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.4);
        hold on;
        plot(timeVals,mean(dataHighStim,1),'Color',[0.7 0.7 0.7],'LineWidth',1.0);
        plot(timeVals,mean(pcaK20Out.cleanedData,1),'LineWidth',1.2);

        tsLegend = {'No-stim','Raw high-stim','PCATemplate K=20'};

        if isstruct(smartaLiteOut) && isfield(smartaLiteOut,'cleanedData')
            plot(timeVals,mean(smartaLiteOut.cleanedData,1),'LineWidth',1.2);
            tsLegend{end+1} = 'SMARTALite';
        end

        if isstruct(smartaFullOut) && isfield(smartaFullOut,'cleanedData')
            plot(timeVals,mean(smartaFullOut.cleanedData,1),'LineWidth',1.2);
            tsLegend{end+1} = 'SMARTAFull';
        end

        xlim([-0.05 0.45]);
        xlabel('Time (s)');
        ylabel('Mean LFP');
        title(['demo16 SMARTA time-series comparison, elec' num2str(elecNum)]);
        legend(tsLegend,'Location','best');
        grid on;

        saveas(figTS,fullfile(figFolder,['demo16_elec' num2str(elecNum) '_time_series.png']));
        savefig(figTS,fullfile(figFolder,['demo16_elec' num2str(elecNum) '_time_series.fig']));
    end

    clear fftSMARTALite fftSMARTAFull
end

%% ========================================================================
% Results table
% ========================================================================

if isempty(resultRows)
    error('No electrodes were processed successfully.');
end

resultsTable = struct2table(resultRows);

fprintf('\nDemo 16 SMARTA results table:\n');
disp(resultsTable);

%% ========================================================================
% Summary statistics
% ========================================================================

errorMatrix = [ ...
    resultsTable.erpError, ...
    resultsTable.alignedError, ...
    resultsTable.hybridError, ...
    resultsTable.pcaK3Error, ...
    resultsTable.pcaK20Error, ...
    resultsTable.smartaLiteError, ...
    resultsTable.smartaFullError];

improvementMatrix = [ ...
    resultsTable.erpImprovement, ...
    resultsTable.alignedImprovement, ...
    resultsTable.hybridImprovement, ...
    resultsTable.pcaK3Improvement, ...
    resultsTable.pcaK20Improvement, ...
    resultsTable.smartaLiteImprovement, ...
    resultsTable.smartaFullImprovement];

meanError = mean(errorMatrix,1,'omitnan');
stdError = std(errorMatrix,0,1,'omitnan');

meanImprovement = mean(improvementMatrix,1,'omitnan');
stdImprovement = std(improvementMatrix,0,1,'omitnan');

summaryTable = table( ...
    methodNames(:), ...
    meanError(:), ...
    stdError(:), ...
    meanImprovement(:), ...
    stdImprovement(:), ...
    'VariableNames',{'method','meanError','stdError','meanImprovement','stdImprovement'} ...
);

summaryTable = sortrows(summaryTable,'meanError','ascend');

fprintf('\nDemo 16 SMARTA summary table:\n');
disp(summaryTable);

%% ========================================================================
% Best method per electrode
% ========================================================================

[~,bestMethodIdx] = min(errorMatrix,[],2,'omitnan');
bestMethodPerElectrode = methodNames(bestMethodIdx)';

bestMethodTable = table( ...
    resultsTable.electrode, ...
    bestMethodPerElectrode, ...
    'VariableNames',{'electrode','bestMethod'} ...
);

fprintf('\nBest method per electrode:\n');
disp(bestMethodTable);

fprintf('\nBest method counts:\n');

for iMethod = 1:length(methodNames)
    nBest = sum(bestMethodIdx == iMethod);
    fprintf('%s: %d electrodes\n',methodNames{iMethod},nBest);
end

%% ========================================================================
% Summary figure
% ========================================================================

figure('Name','demo16 SMARTA summary: improvement');

bar(meanImprovement);
hold on;
errorbar(1:length(methodNames),meanImprovement,stdImprovement,'k.','LineWidth',1.2);

xticks(1:length(methodNames));
xticklabels(methodNames);
xtickangle(30);

ylabel('Improvement relative to raw high-stim (%)');
title('demo16 SMARTA comparison on selected electrodes');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo16_SMARTA_summary_improvement.png'));
    savefig(gcf,fullfile(figFolder,'demo16_SMARTA_summary_improvement.fig'));
end

%% ========================================================================
% Save results
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo16_SMARTA_results_' timestampString '.mat']);
resultsCsvFile = fullfile(resultsFolder,['demo16_SMARTA_results_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo16_SMARTA_summary_' timestampString '.csv']);

save(resultsMatFile, ...
    'resultsTable', ...
    'summaryTable', ...
    'bestMethodTable', ...
    'methodNames', ...
    'electrodesToTest', ...
    'pulseTimes', ...
    'smartaLiteOpts', ...
    'smartaFullOpts');

writetable(resultsTable,resultsCsvFile);
writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved demo_16 results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',resultsCsvFile);
fprintf('%s\n',summaryCsvFile);

if saveFigures
    fprintf('\nSaved demo_16 figures in:\n');
    fprintf('%s\n',figFolder);
end

%% ========================================================================
% Final note
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_16 SMARTA complete\n');
fprintf('============================================================\n');
fprintf('This is a SMARTA-inspired pulse-wise ICMS adaptation.\n');
fprintf('Compare SMARTALite/SMARTAFull against PCATemplate K=20.\n');
fprintf('If SMARTA improves the two exception electrodes elec8 and elec24,\n');
fprintf('it may be worth developing further.\n');