clear; close all; clc;

%% demo_15_RepresentativeElectrodeFigures
%
% Goal:

%
% demo_14 showed that PCATemplate K=20 was best on 29/31 good V1 electrodes.
% This demo plots representative electrodes:
%   - electrodes where PCA K=20 performs strongly
%   - exceptions from demo_14: elec8 and elec24
%
% Figures:
% 1. Mean time-series comparison
% 2. Zoomed artifact-window comparison
% 3. FFT comparison against no-stim reference
% 4. Error bar summary for selected electrodes

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

% Representative electrodes
% elec3,6,19,35: PCA K=20 strong examples
% elec8: TunedHybrid was best in demo_14
% elec24: ERPSubtraction was best in demo_14
electrodesToPlot = [3 6 19 35 8 24];

% Windows
artifactWindow = [-0.02 0.4];
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% PCA settings
pcaKConservative = 3;
pcaKAggressive = 20;

% Save figures?
saveFigures = true;

%% ========================================================================
% Paths and data loading
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
% Confirm good V1 electrode list from RFData
% ========================================================================

v1Electrodes = 1:48;
stimElectrode = 1;

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
goodV1Electrodes = intersect(highRMSElectrodes,v1Electrodes);
goodV1Electrodes = setdiff(goodV1Electrodes,stimElectrode);

fprintf('\nGood V1 electrodes available after excluding elec1:\n');
disp(goodV1Electrodes);

missingElectrodes = setdiff(electrodesToPlot,goodV1Electrodes);
if ~isempty(missingElectrodes)
    warning('These requested electrodes are not in the good V1 list:');
    disp(missingElectrodes);
end

electrodesToPlot = intersect(electrodesToPlot,goodV1Electrodes,'stable');

fprintf('\nElectrodes selected for representative figures:\n');
disp(electrodesToPlot);

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

% Tuned hybrid from demo_14
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
pcaParamsK3.numComponents = pcaKConservative;
pcaParamsK3.removeMeanTemplate = true;
pcaParamsK3.taperEdgeMS = 2;
pcaParamsK3.doBaselineCorrection = false;

% PCA K=20
pcaParamsK20 = struct();
pcaParamsK20.artifactWindow = artifactWindow;
pcaParamsK20.numComponents = pcaKAggressive;
pcaParamsK20.removeMeanTemplate = true;
pcaParamsK20.taperEdgeMS = 2;
pcaParamsK20.doBaselineCorrection = false;

methodNames = { ...
    'Raw high-stim', ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'TunedHybrid', ...
    'PCATemplate K=3', ...
    'PCATemplate K=20'};

%% ========================================================================
% Output folder
% ========================================================================

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo15');

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Loop through representative electrodes
% ========================================================================

summaryRows = struct([]);
rowCounter = 0;

for iElec = 1:length(electrodesToPlot)

    elecNum = electrodesToPlot(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    fprintf('\n============================================================\n');
    fprintf('Processing representative electrode elec%d\n',elecNum);
    fprintf('============================================================\n');

    if ~exist(elecFile,'file')
        warning('File not found: %s',elecFile);
        continue;
    end

    D = load(elecFile);

    if ~isfield(D,'analogData')
        warning('analogData not found in elec%d.',elecNum);
        continue;
    end

    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    %% Apply methods
    erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);
    alignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);
    hybridOut = ERPAlignedPulsewise(dataHighStim,timeVals,hybridParams);
    pcaK3Out = PCATemplate(dataHighStim,timeVals,pcaParamsK3);
    pcaK20Out = PCATemplate(dataHighStim,timeVals,pcaParamsK20);

    %% FFT summaries
    fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
    fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);
    fftERP = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);
    fftAligned = compute_fft_summary(alignedOut.cleanedData,timeVals,fftIdx);
    fftHybrid = compute_fft_summary(hybridOut.cleanedData,timeVals,fftIdx);
    fftPCAK3 = compute_fft_summary(pcaK3Out.cleanedData,timeVals,fftIdx);
    fftPCAK20 = compute_fft_summary(pcaK20Out.cleanedData,timeVals,fftIdx);

    freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
               fftNoStim.freqAxis <= freqRangeForMetric(2);

    rawError = norm(fftHighRaw.logMeanMagnitude(freqMask) - fftNoStim.logMeanMagnitude(freqMask));
    erpError = norm(fftERP.logMeanMagnitude(freqMask) - fftNoStim.logMeanMagnitude(freqMask));
    alignedError = norm(fftAligned.logMeanMagnitude(freqMask) - fftNoStim.logMeanMagnitude(freqMask));
    hybridError = norm(fftHybrid.logMeanMagnitude(freqMask) - fftNoStim.logMeanMagnitude(freqMask));
    pcaK3Error = norm(fftPCAK3.logMeanMagnitude(freqMask) - fftNoStim.logMeanMagnitude(freqMask));
    pcaK20Error = norm(fftPCAK20.logMeanMagnitude(freqMask) - fftNoStim.logMeanMagnitude(freqMask));

    errorVals = [rawError erpError alignedError hybridError pcaK3Error pcaK20Error];

    improvementVals = 100 * (rawError - errorVals) / rawError;

    [bestError,bestIdx] = min(errorVals(2:end));
    bestIdx = bestIdx + 1;
    bestMethod = methodNames{bestIdx};

    fprintf('Raw error        = %.4f\n',rawError);
    fprintf('ERPSubtraction   = %.4f, improvement %.2f %%\n',erpError,improvementVals(2));
    fprintf('ERPAligned       = %.4f, improvement %.2f %%\n',alignedError,improvementVals(3));
    fprintf('TunedHybrid      = %.4f, improvement %.2f %%\n',hybridError,improvementVals(4));
    fprintf('PCATemplate K=3  = %.4f, improvement %.2f %%\n',pcaK3Error,improvementVals(5));
    fprintf('PCATemplate K=20 = %.4f, improvement %.2f %%\n',pcaK20Error,improvementVals(6));
    fprintf('Best method      = %s, error %.4f\n',bestMethod,bestError);

    %% Store summary
    rowCounter = rowCounter + 1;
    summaryRows(rowCounter).electrode = elecNum;
    summaryRows(rowCounter).rawError = rawError;
    summaryRows(rowCounter).erpError = erpError;
    summaryRows(rowCounter).alignedError = alignedError;
    summaryRows(rowCounter).hybridError = hybridError;
    summaryRows(rowCounter).pcaK3Error = pcaK3Error;
    summaryRows(rowCounter).pcaK20Error = pcaK20Error;
    summaryRows(rowCounter).erpImprovement = improvementVals(2);
    summaryRows(rowCounter).alignedImprovement = improvementVals(3);
    summaryRows(rowCounter).hybridImprovement = improvementVals(4);
    summaryRows(rowCounter).pcaK3Improvement = improvementVals(5);
    summaryRows(rowCounter).pcaK20Improvement = improvementVals(6);
    summaryRows(rowCounter).bestMethod = string(bestMethod);

    %% ====================================================================
    % Figure A: time-series comparison
    % ====================================================================

    figA = figure('Name',['demo15 elec' num2str(elecNum) ' time series']);

    subplot(3,1,1);
    plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.2);
    hold on;
    plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.0);
    xlim([-0.2 0.8]);
    title(['elec' num2str(elecNum) ': No-stim vs raw high-stim']);
    xlabel('Time (s)');
    ylabel('Mean LFP');
    legend('No-stim reference','Raw high-stim','Location','best');
    grid on;

    subplot(3,1,2);
    plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.0);
    hold on;
    plot(timeVals,mean(erpOut.cleanedData,1),'b','LineWidth',1.0);
    plot(timeVals,mean(hybridOut.cleanedData,1),'m','LineWidth',1.0);
    plot(timeVals,mean(pcaK20Out.cleanedData,1),'g','LineWidth',1.2);
    xlim([-0.05 0.5]);
    title('Artifact-window mean signal comparison');
    xlabel('Time (s)');
    ylabel('Mean LFP');
    legend('Raw high-stim','ERPSubtraction','TunedHybrid','PCATemplate K=20','Location','best');
    grid on;

    subplot(3,1,3);
    plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.2);
    hold on;
    plot(timeVals,mean(pcaK3Out.cleanedData,1),'c','LineWidth',1.0);
    plot(timeVals,mean(pcaK20Out.cleanedData,1),'g','LineWidth',1.2);
    xlim([-0.2 0.8]);
    title('PCA cleaned mean vs no-stim reference');
    xlabel('Time (s)');
    ylabel('Mean LFP');
    legend('No-stim reference','PCATemplate K=3','PCATemplate K=20','Location','best');
    grid on;

    if saveFigures
        saveas(figA,fullfile(figFolder,['demo15_elec' num2str(elecNum) '_time_series.png']));
        savefig(figA,fullfile(figFolder,['demo15_elec' num2str(elecNum) '_time_series.fig']));
    end

    %% ====================================================================
    % Figure B: FFT comparison
    % ====================================================================

    figB = figure('Name',['demo15 elec' num2str(elecNum) ' FFT']);

    plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.4);
    hold on;
    plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'r','LineWidth',1.0);
    plot(fftERP.freqAxis,fftERP.logMeanMagnitude,'b','LineWidth',1.0);
    plot(fftHybrid.freqAxis,fftHybrid.logMeanMagnitude,'m','LineWidth',1.0);
    plot(fftPCAK3.freqAxis,fftPCAK3.logMeanMagnitude,'c','LineWidth',1.0);
    plot(fftPCAK20.freqAxis,fftPCAK20.logMeanMagnitude,'g','LineWidth',1.3);

    xlim(freqRangeForMetric);
    title(['elec' num2str(elecNum) ': FFT comparison, best = ' bestMethod]);
    xlabel('Frequency (Hz)');
    ylabel('log10 mean FFT magnitude');
    legend( ...
        'No-stim reference', ...
        'Raw high-stim', ...
        'ERPSubtraction', ...
        'TunedHybrid', ...
        'PCATemplate K=3', ...
        'PCATemplate K=20', ...
        'Location','best');
    grid on;

    if saveFigures
        saveas(figB,fullfile(figFolder,['demo15_elec' num2str(elecNum) '_fft.png']));
        savefig(figB,fullfile(figFolder,['demo15_elec' num2str(elecNum) '_fft.fig']));
    end

    %% ====================================================================
    % Figure C: error bar plot for this electrode
    % ====================================================================

    figC = figure('Name',['demo15 elec' num2str(elecNum) ' errors']);

    bar(errorVals);
    xticks(1:length(methodNames));
    xticklabels(methodNames);
    xtickangle(30);
    ylabel('FFT error vs no-stim');
    title(['elec' num2str(elecNum) ': method error comparison']);
    grid on;

    if saveFigures
        saveas(figC,fullfile(figFolder,['demo15_elec' num2str(elecNum) '_errors.png']));
        savefig(figC,fullfile(figFolder,['demo15_elec' num2str(elecNum) '_errors.fig']));
    end
end

%% ========================================================================
% Summary table and combined figure
% ========================================================================

summaryTable = struct2table(summaryRows);

fprintf('\nRepresentative-electrode summary table:\n');
disp(summaryTable);

figSummary = figure('Name','demo15 representative electrode summary');

barData = [ ...
    summaryTable.erpImprovement, ...
    summaryTable.alignedImprovement, ...
    summaryTable.hybridImprovement, ...
    summaryTable.pcaK3Improvement, ...
    summaryTable.pcaK20Improvement];

bar(categorical(string(summaryTable.electrode)),barData);

xlabel('Electrode');
ylabel('Improvement relative to raw high-stim (%)');
title('Representative electrodes: improvement by method');
legend(methodNames(2:end),'Location','best');
grid on;

if saveFigures
    saveas(figSummary,fullfile(figFolder,'demo15_representative_summary.png'));
    savefig(figSummary,fullfile(figFolder,'demo15_representative_summary.fig'));
end

%% ========================================================================
% Save summary
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

timestampString = datestr(now,'yyyymmdd_HHMMSS');

summaryMatFile = fullfile(resultsFolder,['demo15_representative_electrodes_' timestampString '.mat']);
summaryCsvFile = fullfile(resultsFolder,['demo15_representative_electrodes_' timestampString '.csv']);

save(summaryMatFile,'summaryTable','electrodesToPlot','methodNames','artifactWindow','fftWindow','freqRangeForMetric');
writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved demo_15 summary:\n');
fprintf('%s\n',summaryMatFile);
fprintf('%s\n',summaryCsvFile);

if saveFigures
    fprintf('\nSaved demo_15 figures in:\n');
    fprintf('%s\n',figFolder);
end

%% ========================================================================
% Final note
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_15 complete\n');
fprintf('============================================================\n');
fprintf('Use these figures to visually support demo_14 results.\n');
fprintf('demo_14 showed PCATemplate K=20 was best on 29/31 good V1 electrodes.\n');
fprintf('demo_15 shows examples and exceptions electrode-by-electrode.\n');