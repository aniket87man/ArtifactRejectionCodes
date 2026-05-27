clear; close all; clc;

%% demo_11_PCATemplateValidation
%
% Goal:
% Validate PCATemplate and check possible over-removal.
%
% This demo tests:
%
% 1. High-stim correction performance for different K values.
% 2. No-stim distortion control:
%    Apply the same PCATemplate idea to no-stim data and measure how much it
%    changes signal that should not contain stimulation artifact.
%
% Important:
% This no-stim test is a conservative/aggressiveness control.
% If high K strongly changes no-stim data, then high K may be risky.

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

% K values to test
componentList = [0 1 2 3 5 10 20];

% FFT settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% Plot settings
numTrialsToShow = 20;

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

%% ========================================================================
% FFT indices and artifact indices
% ========================================================================

fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));
artifactIdx = find(timeVals >= artifactWindow(1) & timeVals < artifactWindow(2));

fprintf('\nFFT window has %d samples.\n',length(fftIdx));
fprintf('Window duration = %.4f s.\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));
fprintf('Artifact window has %d samples.\n',length(artifactIdx));
fprintf('Artifact window duration = %.4f s.\n',timeVals(artifactIdx(end))-timeVals(artifactIdx(1)));

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
% FFT summaries for reference methods
% ========================================================================

fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);
fftERPClean = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);
fftAlignedClean = compute_fft_summary(erpAlignedOut.cleanedData,timeVals,fftIdx);
fftHybridClean = compute_fft_summary(hybridOut.cleanedData,timeVals,fftIdx);

freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
           fftNoStim.freqAxis <= freqRangeForMetric(2);

rawHighError = norm(fftHighRaw.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

erpError = norm(fftERPClean.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

alignedError = norm(fftAlignedClean.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

hybridError = norm(fftHybridClean.logMeanMagnitude(freqMask) - ...
                   fftNoStim.logMeanMagnitude(freqMask));

erpImprovement = 100 * (rawHighError - erpError) / rawHighError;
alignedImprovement = 100 * (rawHighError - alignedError) / rawHighError;
hybridImprovement = 100 * (rawHighError - hybridError) / rawHighError;

fprintf('\nReference method results:\n');
fprintf('Raw high-stim FFT error vs no-stim = %.4f\n',rawHighError);
fprintf('ERPSubtraction error               = %.4f, improvement = %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned error                   = %.4f, improvement = %.2f %%\n',alignedError,alignedImprovement);
fprintf('Tuned hybrid error                 = %.4f, improvement = %.2f %%\n',hybridError,hybridImprovement);

%% ========================================================================
% PCATemplate validation sweep
% ========================================================================

results = struct([]);

fprintf('\nRunning PCATemplate validation sweep...\n');

for iK = 1:length(componentList)

    K = componentList(iK);

    pcaParams = struct();
    pcaParams.artifactWindow = artifactWindow;
    pcaParams.numComponents = K;
    pcaParams.removeMeanTemplate = true;
    pcaParams.taperEdgeMS = 2;
    pcaParams.doBaselineCorrection = false;

    %% Apply PCATemplate to high-stim data
    pcaHighOut = PCATemplate(dataHighStim,timeVals,pcaParams);
    fftPCAHigh = compute_fft_summary(pcaHighOut.cleanedData,timeVals,fftIdx);

    pcaHighError = norm(fftPCAHigh.logMeanMagnitude(freqMask) - ...
                        fftNoStim.logMeanMagnitude(freqMask));

    pcaHighImprovement = 100 * (rawHighError - pcaHighError) / rawHighError;

    %% Apply PCATemplate to no-stim data as a distortion/aggressiveness control
    pcaNoStimOut = PCATemplate(dataNoStim,timeVals,pcaParams);
    fftNoStimClean = compute_fft_summary(pcaNoStimOut.cleanedData,timeVals,fftIdx);

    % FFT distortion: how much the no-stim spectrum changes after PCA.
    noStimFFTDistortion = norm(fftNoStimClean.logMeanMagnitude(freqMask) - ...
                               fftNoStim.logMeanMagnitude(freqMask));

    % RMS distortion inside artifact window.
    noStimChange = pcaNoStimOut.cleanedData(:,artifactIdx) - dataNoStim(:,artifactIdx);

    noStimRMSChange = sqrt(mean(noStimChange(:).^2));

    tmpNoStim = dataNoStim(:,artifactIdx);
    noStimOriginalRMS = sqrt(mean(tmpNoStim(:).^2));

    noStimRMSChangePercent = 100 * noStimRMSChange / noStimOriginalRMS;

    % High-stim RMS residual after correction inside artifact window.
    tmpHighStimOriginal = dataHighStim(:,artifactIdx);
    tmpHighStimClean = pcaHighOut.cleanedData(:,artifactIdx);

    highStimOriginalRMS = sqrt(mean(tmpHighStimOriginal(:).^2));
    highStimCleanRMS = sqrt(mean(tmpHighStimClean(:).^2));

    highStimRMSReductionPercent = 100 * ...
    (highStimOriginalRMS - highStimCleanRMS) / highStimOriginalRMS;
    %% Store
    results(iK).K = K;
    results(iK).pcaHighError = pcaHighError;
    results(iK).pcaHighImprovement = pcaHighImprovement;
    results(iK).noStimFFTDistortion = noStimFFTDistortion;
    results(iK).noStimRMSChange = noStimRMSChange;
    results(iK).noStimRMSChangePercent = noStimRMSChangePercent;
    results(iK).highStimRMSReductionPercent = highStimRMSReductionPercent;
    results(iK).pcaHighOut = pcaHighOut;
    results(iK).pcaNoStimOut = pcaNoStimOut;
    results(iK).fftPCAHigh = fftPCAHigh;
    results(iK).fftNoStimClean = fftNoStimClean;

    fprintf('K = %2d | high error %.4f, improvement %.2f %% | no-stim FFT distortion %.4f | no-stim RMS change %.2f %%\n', ...
        K,pcaHighError,pcaHighImprovement,noStimFFTDistortion,noStimRMSChangePercent);

end

%% ========================================================================
% Convert results to table
% ========================================================================

Kcol = zeros(length(results),1);
highErrorCol = zeros(length(results),1);
highImprovementCol = zeros(length(results),1);
noStimFFTDistortionCol = zeros(length(results),1);
noStimRMSChangePercentCol = zeros(length(results),1);
highStimRMSReductionPercentCol = zeros(length(results),1);

for i = 1:length(results)
    Kcol(i) = results(i).K;
    highErrorCol(i) = results(i).pcaHighError;
    highImprovementCol(i) = results(i).pcaHighImprovement;
    noStimFFTDistortionCol(i) = results(i).noStimFFTDistortion;
    noStimRMSChangePercentCol(i) = results(i).noStimRMSChangePercent;
    highStimRMSReductionPercentCol(i) = results(i).highStimRMSReductionPercent;
end

validationTable = table( ...
    Kcol, ...
    highErrorCol, ...
    highImprovementCol, ...
    noStimFFTDistortionCol, ...
    noStimRMSChangePercentCol, ...
    highStimRMSReductionPercentCol, ...
    'VariableNames', { ...
        'K', ...
        'highStimError', ...
        'highStimImprovement', ...
        'noStimFFTDistortion', ...
        'noStimRMSChangePercent', ...
        'highStimRMSReductionPercent' ...
    });

fprintf('\nPCATemplate validation table:\n');
disp(validationTable);

%% ========================================================================
% Choose conservative and aggressive PCA examples for visualization
% ========================================================================

% Conservative candidate: K = 3
conservativeK = 3;

% Aggressive candidate: best high-stim FFT error
[~,bestIdx] = min(highErrorCol);
aggressiveK = Kcol(bestIdx);

conservativeIdx = find(Kcol == conservativeK,1);
aggressiveIdx = bestIdx;

pcaConservativeOut = results(conservativeIdx).pcaHighOut;
pcaAggressiveOut = results(aggressiveIdx).pcaHighOut;

fftConservative = results(conservativeIdx).fftPCAHigh;
fftAggressive = results(aggressiveIdx).fftPCAHigh;

fprintf('\nCandidate methods:\n');
fprintf('Conservative PCA K = %d: error %.4f, improvement %.2f %%, no-stim RMS change %.2f %%\n', ...
    conservativeK, ...
    results(conservativeIdx).pcaHighError, ...
    results(conservativeIdx).pcaHighImprovement, ...
    results(conservativeIdx).noStimRMSChangePercent);

fprintf('Aggressive PCA K = %d: error %.4f, improvement %.2f %%, no-stim RMS change %.2f %%\n', ...
    aggressiveK, ...
    results(aggressiveIdx).pcaHighError, ...
    results(aggressiveIdx).pcaHighImprovement, ...
    results(aggressiveIdx).noStimRMSChangePercent);

%% ========================================================================
% Figure 1: High-stim performance vs K
% ========================================================================

figure('Name','PCATemplate validation: high-stim performance');

subplot(2,1,1);
plot(Kcol,highErrorCol,'o-','LineWidth',1.2);
hold on;
yline(erpError,'b--','LineWidth',1.2);
yline(alignedError,'g--','LineWidth',1.2);
yline(hybridError,'m--','LineWidth',1.2);
xlabel('Number of PCA components, K');
ylabel('FFT error vs no-stim');
title('High-stim correction performance');
legend('PCATemplate','ERPSubtraction','ERPAligned','Tuned hybrid','Location','best');
grid on;

subplot(2,1,2);
plot(Kcol,highImprovementCol,'o-','LineWidth',1.2);
hold on;
yline(erpImprovement,'b--','LineWidth',1.2);
yline(alignedImprovement,'g--','LineWidth',1.2);
yline(hybridImprovement,'m--','LineWidth',1.2);
xlabel('Number of PCA components, K');
ylabel('Improvement (%)');
title('High-stim improvement');
legend('PCATemplate','ERPSubtraction','ERPAligned','Tuned hybrid','Location','best');
grid on;

%% ========================================================================
% Figure 2: No-stim distortion vs K
% ========================================================================

figure('Name','PCATemplate validation: no-stim distortion');

subplot(2,1,1);
plot(Kcol,noStimFFTDistortionCol,'o-','LineWidth',1.2);
xlabel('Number of PCA components, K');
ylabel('No-stim FFT distortion');
title('How much PCATemplate changes no-stim spectrum');
grid on;

subplot(2,1,2);
plot(Kcol,noStimRMSChangePercentCol,'o-','LineWidth',1.2);
xlabel('Number of PCA components, K');
ylabel('No-stim RMS change (%)');
title('How much PCATemplate changes no-stim time-domain signal');
grid on;

%% ========================================================================
% Figure 3: Tradeoff plot
% ========================================================================

figure('Name','PCATemplate validation: performance vs distortion tradeoff');

scatter(noStimRMSChangePercentCol,highImprovementCol,80,Kcol,'filled');
colorbar;
xlabel('No-stim RMS change (%)');
ylabel('High-stim improvement (%)');
title('Tradeoff: high-stim artifact removal vs no-stim distortion');
grid on;

for i = 1:length(Kcol)
    text(noStimRMSChangePercentCol(i),highImprovementCol(i),['  K=' num2str(Kcol(i))]);
end

%% ========================================================================
% Figure 4: FFT comparison conservative vs aggressive K
% ========================================================================

figure('Name','PCATemplate validation: FFT conservative vs aggressive');

plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.3);
hold on;
plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'r','LineWidth',1.2);
plot(fftERPClean.freqAxis,fftERPClean.logMeanMagnitude,'b','LineWidth',1.2);
plot(fftHybridClean.freqAxis,fftHybridClean.logMeanMagnitude,'m','LineWidth',1.2);
plot(fftConservative.freqAxis,fftConservative.logMeanMagnitude,'g','LineWidth',1.3);
plot(fftAggressive.freqAxis,fftAggressive.logMeanMagnitude,'c','LineWidth',1.3);

xlim([0 200]);
title('FFT comparison: conservative vs aggressive PCA');
xlabel('Frequency (Hz)');
ylabel('log10 mean FFT magnitude');

legend( ...
    'No-stim reference', ...
    'High-stim raw', ...
    'ERPSubtraction', ...
    'Tuned hybrid', ...
    ['PCATemplate K=' num2str(conservativeK)], ...
    ['PCATemplate K=' num2str(aggressiveK)], ...
    'Location','best');

%% ========================================================================
% Figure 5: Time-domain comparison conservative vs aggressive K
% ========================================================================

figure('Name','PCATemplate validation: time-domain conservative vs aggressive');

subplot(4,1,1);
plot(timeVals,dataHighStim(trialSubset,:)');
xlim([-0.1 0.5]);
title('Raw high-stim trials');
xlabel('Time (s)');
ylabel('LFP');

subplot(4,1,2);
plot(timeVals,pcaConservativeOut.cleanedData(trialSubset,:)');
xlim([-0.1 0.5]);
title(['PCATemplate cleaned trials, conservative K = ' num2str(conservativeK)]);
xlabel('Time (s)');
ylabel('Cleaned LFP');

subplot(4,1,3);
plot(timeVals,pcaAggressiveOut.cleanedData(trialSubset,:)');
xlim([-0.1 0.5]);
title(['PCATemplate cleaned trials, aggressive K = ' num2str(aggressiveK)]);
xlabel('Time (s)');
ylabel('Cleaned LFP');

subplot(4,1,4);
plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.1);
hold on;
plot(timeVals,mean(pcaConservativeOut.cleanedData,1),'g','LineWidth',1.1);
plot(timeVals,mean(pcaAggressiveOut.cleanedData,1),'c','LineWidth',1.1);
xlim([-0.1 0.5]);
title('Mean comparison');
xlabel('Time (s)');
ylabel('Mean LFP');
legend('Raw high-stim',['K=' num2str(conservativeK)],['K=' num2str(aggressiveK)],'Location','best');

%% ========================================================================
% Figure 6: No-stim distortion examples
% ========================================================================

pcaNoStimConservative = results(conservativeIdx).pcaNoStimOut;
pcaNoStimAggressive = results(aggressiveIdx).pcaNoStimOut;

noStimTrialSubset = 1:min(numTrialsToShow,size(dataNoStim,1));

figure('Name','PCATemplate validation: no-stim distortion examples');

subplot(3,1,1);
plot(timeVals,dataNoStim(noStimTrialSubset,:)');
xlim([-0.1 0.5]);
title('Original no-stim trials');
xlabel('Time (s)');
ylabel('LFP');

subplot(3,1,2);
plot(timeVals,pcaNoStimConservative.cleanedData(noStimTrialSubset,:)');
xlim([-0.1 0.5]);
title(['No-stim after PCATemplate, conservative K = ' num2str(conservativeK)]);
xlabel('Time (s)');
ylabel('Cleaned no-stim');

subplot(3,1,3);
plot(timeVals,pcaNoStimAggressive.cleanedData(noStimTrialSubset,:)');
xlim([-0.1 0.5]);
title(['No-stim after PCATemplate, aggressive K = ' num2str(aggressiveK)]);
xlabel('Time (s)');
ylabel('Cleaned no-stim');

%% ========================================================================
% Final printed interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('PCATemplate validation summary\n');
fprintf('============================================================\n');

fprintf('Artifact window = [%.3f %.3f] s\n',artifactWindow(1),artifactWindow(2));
fprintf('FFT metric range = [%d %d] Hz\n',freqRangeForMetric(1),freqRangeForMetric(2));

fprintf('\nReference methods:\n');
fprintf('ERPSubtraction: error %.4f, improvement %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned: error %.4f, improvement %.2f %%\n',alignedError,alignedImprovement);
fprintf('Tuned hybrid: error %.4f, improvement %.2f %%\n',hybridError,hybridImprovement);

fprintf('\nPCATemplate conservative candidate:\n');
fprintf('K = %d\n',conservativeK);
fprintf('High-stim error = %.4f\n',results(conservativeIdx).pcaHighError);
fprintf('High-stim improvement = %.2f %%\n',results(conservativeIdx).pcaHighImprovement);
fprintf('No-stim FFT distortion = %.4f\n',results(conservativeIdx).noStimFFTDistortion);
fprintf('No-stim RMS change = %.2f %%\n',results(conservativeIdx).noStimRMSChangePercent);

fprintf('\nPCATemplate aggressive candidate:\n');
fprintf('K = %d\n',aggressiveK);
fprintf('High-stim error = %.4f\n',results(aggressiveIdx).pcaHighError);
fprintf('High-stim improvement = %.2f %%\n',results(aggressiveIdx).pcaHighImprovement);
fprintf('No-stim FFT distortion = %.4f\n',results(aggressiveIdx).noStimFFTDistortion);
fprintf('No-stim RMS change = %.2f %%\n',results(aggressiveIdx).noStimRMSChangePercent);

fprintf('\nInterpretation:\n');
fprintf('If higher K improves high-stim error but also strongly changes no-stim data,\n');
fprintf('then high K may be over-removing neural signal.\n');
fprintf('A defensible final K should balance high-stim artifact reduction and no-stim preservation.\n');