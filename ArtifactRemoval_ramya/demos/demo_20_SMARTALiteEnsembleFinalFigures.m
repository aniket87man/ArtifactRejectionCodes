clear; close all; clc;

%% demo_20_SMARTALiteEnsembleFinalFigures
%
% Goal:

% SMARTALite Ensemble K3K5.
%
% Representative electrodes:
%   4, 8, 24, 34 : SMARTALite ensemble strong examples
%   6, 35        : PCATemplate K=20 exception examples
%
% Figures per electrode:
%   1. Mean time-series comparison
%   2. FFT comparison
%   3. Method error bar plot
%
% Combined summary:
%   Representative-electrode improvement comparison

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
electrodesToPlot = [4 8 24 34 6 35];

% FFT/artifact settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];
artifactWindow = [-0.02 0.4];

% SMARTALite full-cycle settings
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

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
% Method parameters
% ========================================================================

% ERPSubtraction
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

% PCATemplate K=20
pcaParamsK20 = struct();
pcaParamsK20.artifactWindow = artifactWindow;
pcaParamsK20.numComponents = 20;
pcaParamsK20.removeMeanTemplate = true;
pcaParamsK20.taperEdgeMS = 2;
pcaParamsK20.doBaselineCorrection = false;

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
    'ERPSubtraction', ...
    'PCATemplate K=20', ...
    'SMARTALite K=3', ...
    'SMARTALite K=5', ...
    'SMARTALite Ensemble K3K5'};

%% ========================================================================
% Output folders
% ========================================================================

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo20_SMARTALiteEnsembleFinal');

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

%% ========================================================================
% Loop over representative electrodes
% ========================================================================

summaryRows = struct([]);
rowCounter = 0;

for iElec = 1:length(electrodesToPlot)

    elecNum = electrodesToPlot(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    fprintf('\n============================================================\n');
    fprintf('Processing elec%d\n',elecNum);
    fprintf('============================================================\n');

    if ~exist(elecFile,'file')
        warning('Could not find %s. Skipping.',elecFile);
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

    %% Run methods

    erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);
    erpCleaned = getCleanedData(erpOut);

    pcaK20Out = PCATemplate(dataHighStim,timeVals,pcaParamsK20);
    pcaK20Cleaned = getCleanedData(pcaK20Out);

    smartaK3Out = SMARTALite(dataHighStim,timeVals,smartaK3Opts);
    smartaK3Cleaned = getCleanedData(smartaK3Out);

    smartaK5Out = SMARTALite(dataHighStim,timeVals,smartaK5Opts);
    smartaK5Cleaned = getCleanedData(smartaK5Out);

    ensembleCleaned = 0.5 * (smartaK3Cleaned + smartaK5Cleaned);

    %% FFT summaries and errors

    fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
    fftRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);
    fftERP = compute_fft_summary(erpCleaned,timeVals,fftIdx);
    fftPCAK20 = compute_fft_summary(pcaK20Cleaned,timeVals,fftIdx);
    fftSMARTAK3 = compute_fft_summary(smartaK3Cleaned,timeVals,fftIdx);
    fftSMARTAK5 = compute_fft_summary(smartaK5Cleaned,timeVals,fftIdx);
    fftEnsemble = compute_fft_summary(ensembleCleaned,timeVals,fftIdx);

    freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
               fftNoStim.freqAxis <= freqRangeForMetric(2);

    rawError = computeMetricErrorFromFFT(fftRaw,fftNoStim,freqMask);
    erpError = computeMetricErrorFromFFT(fftERP,fftNoStim,freqMask);
    pcaK20Error = computeMetricErrorFromFFT(fftPCAK20,fftNoStim,freqMask);
    smartaK3Error = computeMetricErrorFromFFT(fftSMARTAK3,fftNoStim,freqMask);
    smartaK5Error = computeMetricErrorFromFFT(fftSMARTAK5,fftNoStim,freqMask);
    ensembleError = computeMetricErrorFromFFT(fftEnsemble,fftNoStim,freqMask);

    errorVals = [erpError pcaK20Error smartaK3Error smartaK5Error ensembleError];
    improvementVals = 100 * (rawError - errorVals) / rawError;

    [bestError,bestIdx] = min(errorVals);
    bestMethod = methodNames{bestIdx};

    fprintf('Raw error                      = %.4f\n',rawError);
    fprintf('ERPSubtraction                 = %.4f, improvement %.2f %%\n',erpError,improvementVals(1));
    fprintf('PCATemplate K=20               = %.4f, improvement %.2f %%\n',pcaK20Error,improvementVals(2));
    fprintf('SMARTALite K=3                 = %.4f, improvement %.2f %%\n',smartaK3Error,improvementVals(3));
    fprintf('SMARTALite K=5                 = %.4f, improvement %.2f %%\n',smartaK5Error,improvementVals(4));
    fprintf('SMARTALite Ensemble K3K5       = %.4f, improvement %.2f %%\n',ensembleError,improvementVals(5));
    fprintf('Best method                    = %s, error %.4f\n',bestMethod,bestError);

    %% Store summary

    rowCounter = rowCounter + 1;

    summaryRows(rowCounter).electrode = elecNum;
    summaryRows(rowCounter).rawError = rawError;

    summaryRows(rowCounter).erpError = erpError;
    summaryRows(rowCounter).erpImprovement = improvementVals(1);

    summaryRows(rowCounter).pcaK20Error = pcaK20Error;
    summaryRows(rowCounter).pcaK20Improvement = improvementVals(2);

    summaryRows(rowCounter).smartaK3Error = smartaK3Error;
    summaryRows(rowCounter).smartaK3Improvement = improvementVals(3);

    summaryRows(rowCounter).smartaK5Error = smartaK5Error;
    summaryRows(rowCounter).smartaK5Improvement = improvementVals(4);

    summaryRows(rowCounter).ensembleError = ensembleError;
    summaryRows(rowCounter).ensembleImprovement = improvementVals(5);

    summaryRows(rowCounter).bestMethod = string(bestMethod);

    %% ====================================================================
    % Figure 1: time-series comparison
    % ====================================================================

    figTS = figure('Name',['demo20 elec' num2str(elecNum) ' time series']);

    tiledlayout(3,1);

    nexttile;
    plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.4);
    hold on;
    plot(timeVals,mean(dataHighStim,1),'Color',[0.7 0.7 0.7],'LineWidth',1.0);
    xlim([-0.1 0.5]);
    xlabel('Time (s)');
    ylabel('Mean LFP');
    title(['elec' num2str(elecNum) ': no-stim vs raw high-stim']);
    legend({'No-stim reference','Raw high-stim'},'Location','best');
    grid on;

    nexttile;
    plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.4);
    hold on;
    plot(timeVals,mean(erpCleaned,1),'LineWidth',1.0);
    plot(timeVals,mean(pcaK20Cleaned,1),'LineWidth',1.0);
    plot(timeVals,mean(ensembleCleaned,1),'LineWidth',1.3);
    xlim([-0.05 0.45]);
    xlabel('Time (s)');
    ylabel('Mean LFP');
    title('Cleaned mean signals');
    legend({'No-stim','ERPSubtraction','PCATemplate K=20','SMARTALite Ensemble K3K5'},'Location','best');
    grid on;

    nexttile;
    plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.4);
    hold on;
    plot(timeVals,mean(smartaK3Cleaned,1),'LineWidth',1.0);
    plot(timeVals,mean(smartaK5Cleaned,1),'LineWidth',1.0);
    plot(timeVals,mean(ensembleCleaned,1),'LineWidth',1.3);
    xlim([-0.05 0.45]);
    xlabel('Time (s)');
    ylabel('Mean LFP');
    title('SMARTALite variants');
    legend({'No-stim','SMARTALite K=3','SMARTALite K=5','Ensemble K3K5'},'Location','best');
    grid on;

    if saveFigures
        saveas(figTS,fullfile(figFolder,['demo20_elec' num2str(elecNum) '_time_series.png']));
        savefig(figTS,fullfile(figFolder,['demo20_elec' num2str(elecNum) '_time_series.fig']));
    end

    %% ====================================================================
    % Figure 2: FFT comparison
    % ====================================================================

    figFFT = figure('Name',['demo20 elec' num2str(elecNum) ' FFT']);

    plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.5);
    hold on;
    plot(fftRaw.freqAxis,fftRaw.logMeanMagnitude,'Color',[0.7 0.7 0.7],'LineWidth',1.0);
    plot(fftERP.freqAxis,fftERP.logMeanMagnitude,'LineWidth',1.0);
    plot(fftPCAK20.freqAxis,fftPCAK20.logMeanMagnitude,'LineWidth',1.0);
    plot(fftEnsemble.freqAxis,fftEnsemble.logMeanMagnitude,'LineWidth',1.4);

    xlim(freqRangeForMetric);
    xlabel('Frequency (Hz)');
    ylabel('log10 mean FFT magnitude');
    title(['elec' num2str(elecNum) ': FFT comparison, best = ' bestMethod]);
    legend({'No-stim','Raw high-stim','ERPSubtraction','PCATemplate K=20','SMARTALite Ensemble K3K5'}, ...
        'Location','best');
    grid on;

    if saveFigures
        saveas(figFFT,fullfile(figFolder,['demo20_elec' num2str(elecNum) '_fft.png']));
        savefig(figFFT,fullfile(figFolder,['demo20_elec' num2str(elecNum) '_fft.fig']));
    end

    %% ====================================================================
    % Figure 3: error bar plot
    % ====================================================================

    figBar = figure('Name',['demo20 elec' num2str(elecNum) ' method errors']);

    bar(errorVals);
    xticks(1:length(methodNames));
    xticklabels(methodNames);
    xtickangle(30);
    ylabel('FFT error vs no-stim');
    title(['elec' num2str(elecNum) ': method error comparison']);
    grid on;

    if saveFigures
        saveas(figBar,fullfile(figFolder,['demo20_elec' num2str(elecNum) '_errors.png']));
        savefig(figBar,fullfile(figFolder,['demo20_elec' num2str(elecNum) '_errors.fig']));
    end

end

%% ========================================================================
% Summary table and combined figure
% ========================================================================

summaryTable = struct2table(summaryRows);

fprintf('\nDemo 20 representative-electrode summary table:\n');
disp(summaryTable);

barData = [ ...
    summaryTable.erpImprovement, ...
    summaryTable.pcaK20Improvement, ...
    summaryTable.smartaK3Improvement, ...
    summaryTable.smartaK5Improvement, ...
    summaryTable.ensembleImprovement];

figSummary = figure('Name','demo20 representative electrode summary');

bar(categorical(string(summaryTable.electrode)),barData);

xlabel('Electrode');
ylabel('Improvement relative to raw high-stim (%)');
title('Representative electrodes: final method comparison');
legend(methodNames,'Location','best');
grid on;

if saveFigures
    saveas(figSummary,fullfile(figFolder,'demo20_representative_summary.png'));
    savefig(figSummary,fullfile(figFolder,'demo20_representative_summary.fig'));
end

%% ========================================================================
% Save summary
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

summaryMatFile = fullfile(resultsFolder,['demo20_SMARTALite_ensemble_final_figures_' timestampString '.mat']);
summaryCsvFile = fullfile(resultsFolder,['demo20_SMARTALite_ensemble_final_figures_' timestampString '.csv']);

save(summaryMatFile, ...
    'summaryTable', ...
    'electrodesToPlot', ...
    'methodNames', ...
    'pulseTimes', ...
    'smartaStimFreq', ...
    'smartaPrePulse', ...
    'artifactWindow', ...
    'fftWindow', ...
    'freqRangeForMetric');

writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved demo_20 summary:\n');
fprintf('%s\n',summaryMatFile);
fprintf('%s\n',summaryCsvFile);

if saveFigures
    fprintf('\nSaved demo_20 figures in:\n');
    fprintf('%s\n',figFolder);
end

fprintf('\n============================================================\n');
fprintf('demo_20 complete\n');
fprintf('============================================================\n');
fprintf('Final method shown: SMARTALite Ensemble K3K5.\n');

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

function err = computeMetricErrorFromFFT(fftCleaned,fftNoStim,freqMask)

    err = norm(fftCleaned.logMeanMagnitude(freqMask) - ...
               fftNoStim.logMeanMagnitude(freqMask));
end