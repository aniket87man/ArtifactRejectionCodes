clear; close all; clc;

%% demo_12_PCATemplateCompareK
%
% Goal:
% Compare PCATemplate results for different K values visually.
%
% This demo helps answer:
% Does high K really remove artifact better, or does it over-clean / flatten
% the signal?
%
% It compares:
% 1. No-stim reference
% 2. Raw high-stim
% 3. PCATemplate K = 3
% 4. PCATemplate K = 5
% 5. PCATemplate K = 10
% 6. PCATemplate K = 20
%
% Important:
% PCA basis vectors are not templates by themselves.
% The actual artifact model is trial-specific:
%
% artifactModel_trial = meanTemplate + score1*PC1 + ... + scoreK*PCK

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName   = 'dona';
gridType      = 'Microelectrode';
expDate       = '290825';
protocolName  = 'GRF_001';
electrodeName = 'elec1';

% Conditions
noStimCondition   = {1,1,1,5,5,4};
highStimCondition = {7,1,1,5,5,4};

% PCA settings
artifactWindow = [-0.02 0.4];
kList = [3 5 10 20];

% FFT settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% Plot settings
numTrialsToShow = 20;
numDetailedTrialsToShow = 5;

%% ========================================================================
% Load data
% ========================================================================

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFile     = fullfile(baseFolder,'segmentedData','lfp',[electrodeName '.mat']);
lfpInfoFile = fullfile(baseFolder,'segmentedData','lfp','lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

D = load(lfpFile);
I = load(lfpInfoFile);
P = load(paramFile);

analogData = D.analogData;
timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

%% ========================================================================
% Select trials
% ========================================================================

noStimTrials = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

dataNoStim = analogData(noStimTrials,:);
dataHighStim = analogData(highStimTrials,:);

fprintf('No-stim trials   : %d\n',size(dataNoStim,1));
fprintf('High-stim trials : %d\n',size(dataHighStim,1));

trialSubset = 1:min(numTrialsToShow,size(dataHighStim,1));
trialSubsetDetailed = 1:min(numDetailedTrialsToShow,size(dataHighStim,1));

%% ========================================================================
% FFT indices
% ========================================================================

fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));

fprintf('\nFFT window has %d samples.\n',length(fftIdx));
fprintf('Window duration = %.4f s.\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% ========================================================================
% Reference methods
% ========================================================================

% ERPSubtraction baseline
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);

% ERPAligned
alignedParams = struct();
alignedParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
alignedParams.alignWindow = [-0.01 0.03];
alignedParams.maxShiftMS = 10;
alignedParams.doBaselineCorrection = false;

erpAlignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);

% Tuned ERPAlignedPulsewise hybrid
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

hybridOut = ERPAlignedPulsewise(dataHighStim,timeVals,hybridParams);

%% ========================================================================
% FFT summaries for reference signals
% ========================================================================

fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);
fftERPClean = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);
fftAlignedClean = compute_fft_summary(erpAlignedOut.cleanedData,timeVals,fftIdx);
fftHybridClean = compute_fft_summary(hybridOut.cleanedData,timeVals,fftIdx);

freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
           fftNoStim.freqAxis <= freqRangeForMetric(2);

rawError = norm(fftHighRaw.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

erpError = norm(fftERPClean.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

alignedError = norm(fftAlignedClean.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

hybridError = norm(fftHybridClean.logMeanMagnitude(freqMask) - ...
                   fftNoStim.logMeanMagnitude(freqMask));

erpImprovement = 100 * (rawError - erpError) / rawError;
alignedImprovement = 100 * (rawError - alignedError) / rawError;
hybridImprovement = 100 * (rawError - hybridError) / rawError;

fprintf('\nReference method results:\n');
fprintf('Raw high-stim FFT error vs no-stim = %.4f\n',rawError);
fprintf('ERPSubtraction error               = %.4f, improvement = %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned error                   = %.4f, improvement = %.2f %%\n',alignedError,alignedImprovement);
fprintf('Tuned hybrid error                 = %.4f, improvement = %.2f %%\n',hybridError,hybridImprovement);

%% ========================================================================
% Run PCATemplate for each K
% ========================================================================

pcaResults = struct([]);

fprintf('\nRunning PCATemplate for selected K values...\n');

for iK = 1:length(kList)

    K = kList(iK);

    pcaParams = struct();
    pcaParams.artifactWindow = artifactWindow;
    pcaParams.numComponents = K;
    pcaParams.removeMeanTemplate = true;
    pcaParams.taperEdgeMS = 2;
    pcaParams.doBaselineCorrection = false;

    pcaOut = PCATemplate(dataHighStim,timeVals,pcaParams);
    fftPCA = compute_fft_summary(pcaOut.cleanedData,timeVals,fftIdx);

    pcaError = norm(fftPCA.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

    pcaImprovement = 100 * (rawError - pcaError) / rawError;

    pcaResults(iK).K = K;
    pcaResults(iK).pcaOut = pcaOut;
    pcaResults(iK).fftPCA = fftPCA;
    pcaResults(iK).pcaError = pcaError;
    pcaResults(iK).pcaImprovement = pcaImprovement;

    fprintf('PCATemplate K = %2d -> error %.4f, improvement %.2f %%\n', ...
        K,pcaError,pcaImprovement);

end

%% ========================================================================
% Summary table
% ========================================================================

Kcol = zeros(length(kList),1);
pcaErrorCol = zeros(length(kList),1);
pcaImprovementCol = zeros(length(kList),1);

for iK = 1:length(kList)
    Kcol(iK) = pcaResults(iK).K;
    pcaErrorCol(iK) = pcaResults(iK).pcaError;
    pcaImprovementCol(iK) = pcaResults(iK).pcaImprovement;
end

pcaCompareTable = table(Kcol,pcaErrorCol,pcaImprovementCol, ...
    'VariableNames',{'K','pcaError','pcaImprovement'});

fprintf('\nPCATemplate K comparison table:\n');
disp(pcaCompareTable);

%% ========================================================================
% Figure 1: Mean time-series comparison
% ========================================================================

figure('Name','PCATemplate Compare K: Mean time series');

numRows = 2 + length(kList);

subplot(numRows,1,1);
plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.2);
xlim([-0.1 0.5]);
title('No-stim reference mean');
xlabel('Time (s)');
ylabel('LFP');

subplot(numRows,1,2);
plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.2);
xlim([-0.1 0.5]);
title('Raw high-stim mean');
xlabel('Time (s)');
ylabel('LFP');

for iK = 1:length(kList)

    subplot(numRows,1,iK+2);
    plot(timeVals,mean(pcaResults(iK).pcaOut.cleanedData,1),'LineWidth',1.2);
    xlim([-0.1 0.5]);
    title(['PCATemplate cleaned mean, K = ' num2str(kList(iK))]);
    xlabel('Time (s)');
    ylabel('LFP');

end

%% ========================================================================
% Figure 2: Trial-level comparison, same selected trials
% ========================================================================

figure('Name','PCATemplate Compare K: Trial-level cleaned data');

numRows = 1 + length(kList);

subplot(numRows,1,1);
plot(timeVals,dataHighStim(trialSubset,:)');
xlim([-0.1 0.5]);
title('Raw high-stim trials');
xlabel('Time (s)');
ylabel('LFP');

for iK = 1:length(kList)

    subplot(numRows,1,iK+1);
    plot(timeVals,pcaResults(iK).pcaOut.cleanedData(trialSubset,:)');
    xlim([-0.1 0.5]);
    title(['Cleaned high-stim trials, PCATemplate K = ' num2str(kList(iK))]);
    xlabel('Time (s)');
    ylabel('Cleaned LFP');

end

%% ========================================================================
% Figure 3: Actual artifact models for same trials
% ========================================================================

figure('Name','PCATemplate Compare K: Artifact models');

for iK = 1:length(kList)

    subplot(length(kList),1,iK);
    plot(timeVals,pcaResults(iK).pcaOut.artifactModel(trialSubset,:)');
    xlim([-0.1 0.5]);
    title(['Trial-specific PCA artifact models, K = ' num2str(kList(iK))]);
    xlabel('Time (s)');
    ylabel('Artifact model');

end

%% ========================================================================
% Figure 4: One detailed trial, raw/model/cleaned for each K
% ========================================================================

detailedTrial = trialSubsetDetailed(1);

figure('Name','PCATemplate Compare K: One detailed trial');

numRows = length(kList);

for iK = 1:length(kList)

    subplot(numRows,1,iK);

    plot(timeVals,dataHighStim(detailedTrial,:),'r','LineWidth',1.0);
    hold on;
    plot(timeVals,pcaResults(iK).pcaOut.artifactModel(detailedTrial,:),'k','LineWidth',1.2);
    plot(timeVals,pcaResults(iK).pcaOut.cleanedData(detailedTrial,:),'b','LineWidth',1.0);

    xlim([-0.1 0.5]);
    title(['Trial ' num2str(detailedTrial) ': raw, artifact model, cleaned, K = ' num2str(kList(iK))]);
    xlabel('Time (s)');
    ylabel('LFP');

end

legend('Raw high-stim','Trial-specific artifact model','Cleaned trial','Location','best');

%% ========================================================================
% Figure 5: Difference between K values
% ========================================================================
%
% This shows how much extra signal is removed when increasing K.

figure('Name','PCATemplate Compare K: Extra removed signal');

baseIdx = find(kList == 3,1);

if isempty(baseIdx)
    baseIdx = 1;
end

baseK = kList(baseIdx);
baseArtifactModel = pcaResults(baseIdx).pcaOut.artifactModel;

for iK = 1:length(kList)

    extraRemoved = pcaResults(iK).pcaOut.artifactModel - baseArtifactModel;

    subplot(length(kList),1,iK);
    plot(timeVals,extraRemoved(trialSubset,:)');
    xlim([-0.1 0.5]);
    title(['Extra removed relative to K = ' num2str(baseK) ', current K = ' num2str(kList(iK))]);
    xlabel('Time (s)');
    ylabel('Extra model');

end

%% ========================================================================
% Figure 6: FFT comparison
% ========================================================================

figure('Name','PCATemplate Compare K: FFT');

plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.4);
hold on;
plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'r','LineWidth',1.2);
plot(fftERPClean.freqAxis,fftERPClean.logMeanMagnitude,'b','LineWidth',1.1);
plot(fftHybridClean.freqAxis,fftHybridClean.logMeanMagnitude,'m','LineWidth',1.1);

legendEntries = { ...
    'No-stim reference', ...
    'High-stim raw', ...
    'ERPSubtraction', ...
    'Tuned hybrid'};

for iK = 1:length(kList)

    plot(pcaResults(iK).fftPCA.freqAxis, ...
         pcaResults(iK).fftPCA.logMeanMagnitude, ...
         'LineWidth',1.2);

    legendEntries{end+1} = ['PCATemplate K=' num2str(kList(iK))]; %#ok<SAGROW>

end

xlim([0 200]);
title('FFT comparison across K values');
xlabel('Frequency (Hz)');
ylabel('log10 mean FFT magnitude');
legend(legendEntries,'Location','best');

%% ========================================================================
% Figure 7: Error and improvement vs K
% ========================================================================

figure('Name','PCATemplate Compare K: Metrics');

subplot(2,1,1);
plot(Kcol,pcaErrorCol,'o-','LineWidth',1.2);
hold on;
yline(erpError,'b--','LineWidth',1.2);
yline(alignedError,'g--','LineWidth',1.2);
yline(hybridError,'m--','LineWidth',1.2);
xlabel('Number of PCA components, K');
ylabel('FFT error vs no-stim');
title('PCATemplate FFT error vs K');
legend('PCATemplate','ERPSubtraction','ERPAligned','Tuned hybrid','Location','best');
grid on;

subplot(2,1,2);
plot(Kcol,pcaImprovementCol,'o-','LineWidth',1.2);
hold on;
yline(erpImprovement,'b--','LineWidth',1.2);
yline(alignedImprovement,'g--','LineWidth',1.2);
yline(hybridImprovement,'m--','LineWidth',1.2);
xlabel('Number of PCA components, K');
ylabel('Improvement (%)');
title('PCATemplate improvement vs K');
legend('PCATemplate','ERPSubtraction','ERPAligned','Tuned hybrid','Location','best');
grid on;

%% ========================================================================
% Figure 8: Explained variance comparison
% ========================================================================

figure('Name','PCATemplate Compare K: Explained variance');

% Use the largest K run, but explained variance is computed from the full SVD.
lastPCAOut = pcaResults(end).pcaOut;

numToShow = min(30,length(lastPCAOut.explainedVariance));

subplot(2,1,1);
bar(1:numToShow,100*lastPCAOut.explainedVariance(1:numToShow));
xlabel('PCA component');
ylabel('Explained variance (%)');
title('Explained variance of PCA components');

subplot(2,1,2);
plot(1:numToShow,100*cumsum(lastPCAOut.explainedVariance(1:numToShow)),'o-','LineWidth',1.2);
xlabel('Number of PCA components');
ylabel('Cumulative explained variance (%)');
title('Cumulative explained variance');
grid on;

%% ========================================================================
% Final printed interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('PCATemplate Compare K summary\n');
fprintf('============================================================\n');

fprintf('Artifact window = [%.3f %.3f] s\n',artifactWindow(1),artifactWindow(2));
fprintf('FFT metric range = [%d %d] Hz\n',freqRangeForMetric(1),freqRangeForMetric(2));

fprintf('\nReference methods:\n');
fprintf('ERPSubtraction: error %.4f, improvement %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned: error %.4f, improvement %.2f %%\n',alignedError,alignedImprovement);
fprintf('Tuned hybrid: error %.4f, improvement %.2f %%\n',hybridError,hybridImprovement);

fprintf('\nPCATemplate candidates:\n');

for iK = 1:length(kList)
    fprintf('K = %2d: error %.4f, improvement %.2f %%\n', ...
        kList(iK),pcaErrorCol(iK),pcaImprovementCol(iK));
end

[bestError,bestIdx] = min(pcaErrorCol);
bestK = kList(bestIdx);
bestImprovement = pcaImprovementCol(bestIdx);

fprintf('\nBest K by FFT metric = %d\n',bestK);
fprintf('Best PCA error = %.4f\n',bestError);
fprintf('Best PCA improvement = %.2f %%\n',bestImprovement);

fprintf('\nInterpretation reminder:\n');
fprintf('K=3 or K=5 is more conservative and easier to justify.\n');
fprintf('K=20 gives best FFT metric but may be over-aggressive.\n');
fprintf('Use trial-level figures to check whether K=20 removes plausible neural activity.\n');