clear; close all; clc;

%% demo_09_PCATemplate
% Compare:
% 1. No-stim reference
% 2. High-stim raw
% 3. ERPSubtraction baseline
% 4. ERPAligned
% 5. Best tuned ERPAlignedPulsewise hybrid
% 6. PCATemplate with different numbers of PCA components

%% Paths
folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName  = 'dona';
gridType     = 'Microelectrode';
expDate      = '290825';
protocolName = 'GRF_001';
electrodeName = 'elec1';

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

%% Load data
lfpFile     = fullfile(baseFolder,'segmentedData','lfp',[electrodeName '.mat']);
lfpInfoFile = fullfile(baseFolder,'segmentedData','lfp','lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

D = load(lfpFile);
I = load(lfpInfoFile);
P = load(paramFile);

analogData = D.analogData;
timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

%% Conditions
noStimTrials = parameterCombinations{1,1,1,5,5,4};
highStimTrials = parameterCombinations{7,1,1,5,5,4};

dataNoStim = analogData(noStimTrials,:);
dataHighStim = analogData(highStimTrials,:);

fprintf('No-stim trials   : %d\n',size(dataNoStim,1));
fprintf('High-stim trials : %d\n',size(dataHighStim,1));

%% FFT window and metric range
fftWindow = [0 0.4];
fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));

freqRangeForMetric = [0 200];

fprintf('\nFFT window has %d samples.\n',length(fftIdx));
fprintf('Window duration = %.4f s.\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% Existing methods for comparison

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

% Best tuned ERPAlignedPulsewise from previous sweep:
% pulseWindowMS = [-10 40], thresholdMAD = 2, templateStatistic = median
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

%% FFT summaries for reference methods
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
fprintf('Raw FFT error vs no-stim          = %.4f\n',rawError);
fprintf('ERPSubtraction error              = %.4f, improvement = %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned error                  = %.4f, improvement = %.2f %%\n',alignedError,alignedImprovement);
fprintf('Tuned ERPAlignedPulsewise error   = %.4f, improvement = %.2f %%\n',hybridError,hybridImprovement);

%% PCA component sweep
componentList = [0 1 2 3 5 10 20];

pcaResults = struct([]);

fprintf('\nRunning PCATemplate component sweep...\n');

for iComp = 1:length(componentList)

    pcaParams = struct();

    % Include the first pulse near -10 ms and the full early artifact period.
    pcaParams.artifactWindow = [-0.02 0.4];

    pcaParams.numComponents = componentList(iComp);
    pcaParams.removeMeanTemplate = true;
    pcaParams.taperEdgeMS = 2;
    pcaParams.doBaselineCorrection = false;

    pcaOut = PCATemplate(dataHighStim,timeVals,pcaParams);

    fftPCAClean = compute_fft_summary(pcaOut.cleanedData,timeVals,fftIdx);

    pcaError = norm(fftPCAClean.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

    pcaImprovement = 100 * (rawError - pcaError) / rawError;

    pcaResults(iComp).numComponents = componentList(iComp);
    pcaResults(iComp).pcaError = pcaError;
    pcaResults(iComp).pcaImprovement = pcaImprovement;
    pcaResults(iComp).pcaOut = pcaOut;
    pcaResults(iComp).fftPCAClean = fftPCAClean;

    fprintf('PCATemplate K = %2d -> error %.4f, improvement %.2f %%\n', ...
        componentList(iComp),pcaError,pcaImprovement);
end

%% Convert PCA results to table
numComponentsCol = zeros(length(pcaResults),1);
pcaErrorCol = zeros(length(pcaResults),1);
pcaImprovementCol = zeros(length(pcaResults),1);

for i = 1:length(pcaResults)
    numComponentsCol(i) = pcaResults(i).numComponents;
    pcaErrorCol(i) = pcaResults(i).pcaError;
    pcaImprovementCol(i) = pcaResults(i).pcaImprovement;
end

pcaTable = table(numComponentsCol,pcaErrorCol,pcaImprovementCol, ...
    'VariableNames',{'numComponents','pcaError','pcaImprovement'});

pcaTable = sortrows(pcaTable,'pcaError','ascend');

fprintf('\nPCATemplate results sorted by error:\n');
disp(pcaTable);

%% Best PCA setting
bestK = pcaTable.numComponents(1);
bestIdx = find(componentList == bestK,1);

bestPCAOut = pcaResults(bestIdx).pcaOut;
fftBestPCAClean = pcaResults(bestIdx).fftPCAClean;

fprintf('\nBest PCATemplate setting:\n');
fprintf('numComponents = %d\n',bestK);
fprintf('PCA error     = %.4f\n',pcaResults(bestIdx).pcaError);
fprintf('PCA improvement = %.2f %%\n',pcaResults(bestIdx).pcaImprovement);

%% Figure 1: PCA error vs number of components
figure;

plot(numComponentsCol,pcaErrorCol,'o-','LineWidth',1.2);
hold on;
yline(erpError,'b--','LineWidth',1.2);
yline(alignedError,'g--','LineWidth',1.2);
yline(hybridError,'m--','LineWidth',1.2);

xlabel('Number of PCA components');
ylabel('FFT error vs no-stim');
title('PCATemplate component sweep');
legend('PCATemplate','ERPSubtraction','ERPAligned','Tuned hybrid','Location','best');
grid on;

%% Figure 2: PCA improvement vs number of components
figure;

plot(numComponentsCol,pcaImprovementCol,'o-','LineWidth',1.2);
hold on;
yline(erpImprovement,'b--','LineWidth',1.2);
yline(alignedImprovement,'g--','LineWidth',1.2);
yline(hybridImprovement,'m--','LineWidth',1.2);

xlabel('Number of PCA components');
ylabel('Improvement relative to raw high-stim (%)');
title('PCATemplate improvement vs number of components');
legend('PCATemplate','ERPSubtraction','ERPAligned','Tuned hybrid','Location','best');
grid on;

%% Figure 3: explained variance from best PCA run
figure;

numToShow = min(20,length(bestPCAOut.explainedVariance));

bar(1:numToShow,100*bestPCAOut.explainedVariance(1:numToShow));
xlabel('PCA component');
ylabel('Explained variance (%)');
title('PCATemplate explained variance');

%% Figure 4: PC templates from best PCA run
figure;

numTemplatesToShow = min(5,size(bestPCAOut.pcTemplates,1));

for i = 1:numTemplatesToShow
    subplot(numTemplatesToShow,1,i);

    artifactTime = timeVals(bestPCAOut.artifactIdx);
    plot(artifactTime,bestPCAOut.pcTemplates(i,:),'LineWidth',1.2);

    xlim([-0.05 0.45]);
    ylabel(['PC ' num2str(i)]);

    if i == 1
        title('Top PCA artifact templates');
    end

    if i == numTemplatesToShow
        xlabel('Time (s)');
    end
end

%% Figure 5: trial-level diagnostic for best PCA
numTrialsToShow = 20;
trialSubset = 1:min(numTrialsToShow,size(dataHighStim,1));

figure;

subplot(4,1,1);
plot(timeVals,dataHighStim(trialSubset,:)');
xlim([-0.2 0.8]);
title('High-stim trials before correction');
xlabel('Time (s)');
ylabel('LFP');

subplot(4,1,2);
plot(timeVals,mean(bestPCAOut.artifactModel(trialSubset,:),1),'k','LineWidth',1.5);
xlim([-0.2 0.8]);
title('Average PCA artifact model for selected trials');
xlabel('Time (s)');
ylabel('Artifact model');

subplot(4,1,3);
plot(timeVals,bestPCAOut.cleanedData(trialSubset,:)');
xlim([-0.2 0.8]);
title('High-stim trials after PCATemplate correction');
xlabel('Time (s)');
ylabel('Cleaned LFP');

subplot(4,1,4);
plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.1);
hold on;
plot(timeVals,mean(bestPCAOut.cleanedData,1),'k','LineWidth',1.1);
xlim([-0.2 0.8]);
title('Mean high-stim signal before and after PCATemplate');
xlabel('Time (s)');
ylabel('Mean LFP');
legend('Raw high-stim','PCATemplate cleaned','Location','best');

%% Figure 6: FFT comparison
figure;

plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.3);
hold on;
plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'r','LineWidth',1.2);
plot(fftERPClean.freqAxis,fftERPClean.logMeanMagnitude,'b','LineWidth',1.2);
plot(fftAlignedClean.freqAxis,fftAlignedClean.logMeanMagnitude,'g','LineWidth',1.3);
plot(fftHybridClean.freqAxis,fftHybridClean.logMeanMagnitude,'m','LineWidth',1.3);
plot(fftBestPCAClean.freqAxis,fftBestPCAClean.logMeanMagnitude,'c','LineWidth',1.3);

xlim([0 200]);
title('FFT comparison: PCATemplate vs previous best methods');
xlabel('Frequency (Hz)');
ylabel('log10 mean FFT magnitude');

legend( ...
    'No-stim reference', ...
    'High-stim raw', ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'Tuned ERPAlignedPulsewise', ...
    ['Best PCATemplate K=' num2str(bestK)], ...
    'Location','best');

%% Final summary
fprintf('\nFinal PCATemplate summary:\n');
fprintf('Raw FFT error              = %.4f\n',rawError);
fprintf('ERPSubtraction error       = %.4f, improvement = %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned error           = %.4f, improvement = %.2f %%\n',alignedError,alignedImprovement);
fprintf('Tuned hybrid error         = %.4f, improvement = %.2f %%\n',hybridError,hybridImprovement);
fprintf('Best PCATemplate error     = %.4f, improvement = %.2f %%\n', ...
    pcaResults(bestIdx).pcaError,pcaResults(bestIdx).pcaImprovement);
fprintf('Best PCATemplate K         = %d\n',bestK);