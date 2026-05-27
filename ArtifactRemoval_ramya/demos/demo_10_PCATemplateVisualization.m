clear; close all; clc;

%% demo_10_PCATemplateVisualization
%
% Goal:
% Visualize exactly what PCATemplate is doing.
%
% Important PCA interpretation:
%
% PCA basis vectors are NOT templates by themselves.
%
% For each trial:
%
% artifactModel_trial = meanTemplate ...
%                     + score1 * PC1 ...
%                     + score2 * PC2 ...
%                     + ... + scoreK * PCK
%
% So the true PCA artifact estimate is trial-specific.

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
selectedK = 3;                % Start with 3 for visualization
artifactWindow = [-0.02 0.4]; % PCA artifact window

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
% Select no-stim and high-stim trials
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
% Run reference methods for comparison
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

% Tuned ERPAlignedPulsewise hybrid from previous sweep
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
% Run PCATemplate for selected K
% ========================================================================

pcaParams = struct();
pcaParams.artifactWindow = artifactWindow;
pcaParams.numComponents = selectedK;
pcaParams.removeMeanTemplate = true;
pcaParams.taperEdgeMS = 2;
pcaParams.doBaselineCorrection = false;

pcaOut = PCATemplate(dataHighStim,timeVals,pcaParams);

K = pcaParams.numComponents;

artifactIdx = pcaOut.artifactIdx;
artifactTime = timeVals(artifactIdx);
meanTemplateWindow = pcaOut.meanTemplate(artifactIdx);

fprintf('\nRunning PCATemplate visualization with K = %d\n',selectedK);
fprintf('Artifact window = [%.4f %.4f] s\n',artifactWindow(1),artifactWindow(2));

%% ========================================================================
% FFT summaries and metrics
% ========================================================================

fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);
fftERPClean = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);
fftAlignedClean = compute_fft_summary(erpAlignedOut.cleanedData,timeVals,fftIdx);
fftHybridClean = compute_fft_summary(hybridOut.cleanedData,timeVals,fftIdx);
fftPCAClean = compute_fft_summary(pcaOut.cleanedData,timeVals,fftIdx);

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

pcaError = norm(fftPCAClean.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

erpImprovement = 100 * (rawError - erpError) / rawError;
alignedImprovement = 100 * (rawError - alignedError) / rawError;
hybridImprovement = 100 * (rawError - hybridError) / rawError;
pcaImprovement = 100 * (rawError - pcaError) / rawError;

fprintf('\nMetric results for selected K:\n');
fprintf('Raw FFT error              = %.4f\n',rawError);
fprintf('ERPSubtraction error       = %.4f, improvement = %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned error           = %.4f, improvement = %.2f %%\n',alignedError,alignedImprovement);
fprintf('Tuned hybrid error         = %.4f, improvement = %.2f %%\n',hybridError,hybridImprovement);
fprintf('PCATemplate K=%d error      = %.4f, improvement = %.2f %%\n',selectedK,pcaError,pcaImprovement);

%% ========================================================================
% Figure 1: Raw trials, true PCA artifact model, cleaned trials
% ========================================================================

figure('Name','PCATemplate: Raw, artifact model, cleaned');

subplot(3,1,1);
plot(timeVals,dataHighStim(trialSubset,:)');
xlim([-0.1 0.5]);
title('Raw high-stim trials');
xlabel('Time (s)');
ylabel('LFP');

subplot(3,1,2);
plot(timeVals,pcaOut.artifactModel(trialSubset,:)');
xlim([-0.1 0.5]);
title(['Trial-specific PCA artifact models actually subtracted, K = ' num2str(K)]);
xlabel('Time (s)');
ylabel('Artifact model');

subplot(3,1,3);
plot(timeVals,pcaOut.cleanedData(trialSubset,:)');
xlim([-0.1 0.5]);
title(['Cleaned high-stim trials after PCATemplate, K = ' num2str(K)]);
xlabel('Time (s)');
ylabel('Cleaned LFP');

%% ========================================================================
% Figure 2: Mean template and PCA basis vectors
% ========================================================================
%
% These PC shapes are basis vectors, not artifact templates by themselves.

figure('Name','PCATemplate: Mean template and PCA basis vectors');

subplot(K+1,1,1);
plot(artifactTime,meanTemplateWindow,'k','LineWidth',1.5);
xlim(artifactWindow);
title('Mean artifact template');
xlabel('Time (s)');
ylabel('Amplitude');

for iPC = 1:K

    subplot(K+1,1,iPC+1);
    plot(artifactTime,pcaOut.pcTemplates(iPC,:),'LineWidth',1.2);
    xlim(artifactWindow);
    title(['PC ' num2str(iPC) ' basis vector, not a template']);
    xlabel('Time (s)');
    ylabel(['PC ' num2str(iPC)]);

end

%% ========================================================================
% Figure 3: Weighted PC contributions for selected trials
% ========================================================================
%
% This shows score * PC, which is the actual PCA contribution.

figure('Name','PCATemplate: Weighted PC contributions');

for iTrialPlot = 1:length(trialSubsetDetailed)

    tr = trialSubsetDetailed(iTrialPlot);

    subplot(length(trialSubsetDetailed),1,iTrialPlot);

    plot(artifactTime,meanTemplateWindow,'k','LineWidth',1.2);
    hold on;

    for iPC = 1:K
        pcContribution = pcaOut.scores(tr,iPC) * pcaOut.pcTemplates(iPC,:);
        plot(artifactTime,pcContribution,'LineWidth',1.0);
    end

    xlim(artifactWindow);

    if iTrialPlot == 1
        title('Mean template and weighted PC contributions');
    end

    ylabel(['Trial ' num2str(tr)]);

    if iTrialPlot == length(trialSubsetDetailed)
        xlabel('Time (s)');
    end

end

legendLabels = cell(1,K+1);
legendLabels{1} = 'Mean template';

for iPC = 1:K
    legendLabels{iPC+1} = ['score' num2str(iPC) ' * PC' num2str(iPC)];
end

legend(legendLabels,'Location','best');

%% ========================================================================
% Figure 4: Actual trial-specific PCA subtraction
% ========================================================================
%
% This is the most important figure for explanation:
% raw trial, trial-specific PCA artifact model, cleaned trial.

figure('Name','PCATemplate: Actual trial-specific subtraction');

for iTrialPlot = 1:length(trialSubsetDetailed)

    tr = trialSubsetDetailed(iTrialPlot);

    subplot(length(trialSubsetDetailed),1,iTrialPlot);

    plot(timeVals,dataHighStim(tr,:),'r','LineWidth',1.0);
    hold on;
    plot(timeVals,pcaOut.artifactModel(tr,:),'k','LineWidth',1.3);
    plot(timeVals,pcaOut.cleanedData(tr,:),'b','LineWidth',1.0);

    xlim([-0.1 0.5]);

    if iTrialPlot == 1
        title(['Actual PCATemplate subtraction, K = ' num2str(K)]);
    end

    ylabel(['Trial ' num2str(tr)]);

    if iTrialPlot == length(trialSubsetDetailed)
        xlabel('Time (s)');
    end

end

legend('Raw high-stim trial','Trial-specific PCA artifact model','Cleaned trial','Location','best');

%% ========================================================================
% Figure 5: Cumulative reconstruction for one trial
% ========================================================================
%
% This shows how the artifact model changes as PCs are added.

tr = trialSubsetDetailed(1);

figure('Name','PCATemplate: Cumulative artifact reconstruction');

cumulativeModel = meanTemplateWindow;

subplot(K+1,1,1);
plot(artifactTime,meanTemplateWindow,'k','LineWidth',1.5);
xlim(artifactWindow);
title(['Trial ' num2str(tr) ': mean template only']);
xlabel('Time (s)');
ylabel('Amplitude');

for iPC = 1:K

    cumulativeModel = cumulativeModel + ...
        pcaOut.scores(tr,iPC) * pcaOut.pcTemplates(iPC,:);

    subplot(K+1,1,iPC+1);
    plot(artifactTime,cumulativeModel,'LineWidth',1.2);
    xlim(artifactWindow);
    title(['Trial ' num2str(tr) ': mean template + PCs 1 to ' num2str(iPC)]);
    xlabel('Time (s)');
    ylabel('Amplitude');

end

%% ========================================================================
% Figure 6: Mean view
% ========================================================================

figure('Name','PCATemplate: Mean view');

plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.2);
hold on;
plot(timeVals,mean(pcaOut.artifactModel,1),'k','LineWidth',1.4);
plot(timeVals,mean(pcaOut.cleanedData,1),'b','LineWidth',1.2);

xlim([-0.1 0.5]);
title(['Mean view of PCATemplate subtraction, K = ' num2str(K)]);
xlabel('Time (s)');
ylabel('LFP');
legend('Raw high-stim mean','Mean of trial-specific PCA artifact models','Cleaned mean','Location','best');

%% ========================================================================
% Figure 7: PCA score distributions
% ========================================================================

figure('Name','PCATemplate: PCA score distributions');

for iPC = 1:K

    subplot(K,1,iPC);
    histogram(pcaOut.scores(:,iPC),30);
    title(['Distribution of PCA scores for PC ' num2str(iPC)]);
    xlabel('Score');
    ylabel('Number of trials');

end

%% ========================================================================
% Figure 8: Explained variance
% ========================================================================

figure('Name','PCATemplate: Explained variance');

numToShow = min(30,length(pcaOut.explainedVariance));

bar(1:numToShow,100*pcaOut.explainedVariance(1:numToShow));
xlabel('PCA component');
ylabel('Explained variance (%)');
title('Explained variance of PCA components');

%% ========================================================================
% Figure 9: Cumulative explained variance
% ========================================================================

figure('Name','PCATemplate: Cumulative explained variance');

cumExplained = cumsum(pcaOut.explainedVariance);

plot(1:numToShow,100*cumExplained(1:numToShow),'o-','LineWidth',1.2);
xlabel('Number of PCA components');
ylabel('Cumulative explained variance (%)');
title('Cumulative explained variance');
grid on;

%% ========================================================================
% Figure 10: FFT comparison
% ========================================================================

figure('Name','PCATemplate: FFT comparison');

plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.3);
hold on;
plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'r','LineWidth',1.2);
plot(fftERPClean.freqAxis,fftERPClean.logMeanMagnitude,'b','LineWidth',1.2);
plot(fftAlignedClean.freqAxis,fftAlignedClean.logMeanMagnitude,'g','LineWidth',1.2);
plot(fftHybridClean.freqAxis,fftHybridClean.logMeanMagnitude,'m','LineWidth',1.2);
plot(fftPCAClean.freqAxis,fftPCAClean.logMeanMagnitude,'c','LineWidth',1.3);

xlim([0 200]);
title(['FFT comparison with PCATemplate K = ' num2str(K)]);
xlabel('Frequency (Hz)');
ylabel('log10 mean FFT magnitude');

legend( ...
    'No-stim reference', ...
    'High-stim raw', ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'Tuned ERPAlignedPulsewise', ...
    ['PCATemplate K=' num2str(K)], ...
    'Location','best');

%% ========================================================================
% Figure 11: Quick K sweep
% ========================================================================

componentList = [0 1 2 3 5 10 20];

pcaErrors = zeros(length(componentList),1);
pcaImprovements = zeros(length(componentList),1);

fprintf('\nQuick PCA K sweep:\n');

for iK = 1:length(componentList)

    tempParams = pcaParams;
    tempParams.numComponents = componentList(iK);

    tempOut = PCATemplate(dataHighStim,timeVals,tempParams);
    tempFFT = compute_fft_summary(tempOut.cleanedData,timeVals,fftIdx);

    tempError = norm(tempFFT.logMeanMagnitude(freqMask) - ...
                     fftNoStim.logMeanMagnitude(freqMask));

    tempImprovement = 100 * (rawError - tempError) / rawError;

    pcaErrors(iK) = tempError;
    pcaImprovements(iK) = tempImprovement;

    fprintf('K = %2d -> error %.4f, improvement %.2f %%\n', ...
        componentList(iK),tempError,tempImprovement);

end

figure('Name','PCATemplate: K sweep');

subplot(2,1,1);
plot(componentList,pcaErrors,'o-','LineWidth',1.2);
hold on;
yline(erpError,'b--','LineWidth',1.2);
yline(alignedError,'g--','LineWidth',1.2);
yline(hybridError,'m--','LineWidth',1.2);
xlabel('Number of PCA components, K');
ylabel('FFT error vs no-stim');
title('PCATemplate error vs K');
legend('PCATemplate','ERPSubtraction','ERPAligned','Tuned hybrid','Location','best');
grid on;

subplot(2,1,2);
plot(componentList,pcaImprovements,'o-','LineWidth',1.2);
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
% Final printed interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('PCATemplate visualization summary\n');
fprintf('============================================================\n');

fprintf('Selected K for visualization = %d\n',selectedK);
fprintf('Artifact window              = [%.3f %.3f] s\n',artifactWindow(1),artifactWindow(2));
fprintf('FFT metric range             = [%d %d] Hz\n',freqRangeForMetric(1),freqRangeForMetric(2));

fprintf('\nSelected K metric:\n');
fprintf('PCATemplate K=%d error        = %.4f\n',selectedK,pcaError);
fprintf('PCATemplate K=%d improvement  = %.2f %%\n',selectedK,pcaImprovement);

fprintf('\nQuick K sweep result:\n');

[bestError,bestIdx] = min(pcaErrors);
bestK = componentList(bestIdx);
bestImprovement = pcaImprovements(bestIdx);

fprintf('Best K from sweep             = %d\n',bestK);
fprintf('Best PCA error                = %.4f\n',bestError);
fprintf('Best PCA improvement          = %.2f %%\n',bestImprovement);

fprintf('\nInterpretation reminder:\n');
fprintf('PC shapes are basis vectors, not templates by themselves.\n');
fprintf('The actual PCA template is trial-specific: mean + weighted PCs.\n');
fprintf('Figure 4 is the most important figure for explaining what is subtracted.\n');
fprintf('Higher K may improve FFT error but may also remove neural signal.\n');