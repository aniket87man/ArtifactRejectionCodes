clear; close all; clc;

%% demo_06_PulsewiseTemplate
% Compare:
% 1. No-stim reference
% 2. High-stim raw
% 3. ERPSubtraction baseline
% 4. ERPAligned, previous best method
% 5. PulsewiseTemplate

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

D = load(lfpFile);       % analogData
I = load(lfpInfoFile);   % timeVals
P = load(paramFile);     % parameterCombinations

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

%% FFT window: 0 to 0.4 s
fftWindow = [0 0.4];
fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));

fprintf('\nFFT window has %d samples.\n',length(fftIdx));
fprintf('Window duration = %.4f s.\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% Method 1: ERPSubtraction baseline
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);

%% Method 2: ERPAligned, previous best method
alignedParams = struct();
alignedParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
alignedParams.alignWindow = [-0.01 0.03];
alignedParams.maxShiftMS = 10;
alignedParams.doBaselineCorrection = false;

erpAlignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);

%% Method 3: PulsewiseTemplate
pulseParams = struct();

% Search range where the stimulation pulse train occurs.
pulseParams.pulseSearchWindow = [-0.02 0.32];

% Local window around each detected pulse.
% This subtracts from 5 ms before pulse to 35 ms after pulse.
pulseParams.pulseWindowMS = [-5 35];

% Minimum gap between detected pulses.
pulseParams.minPulseDistanceMS = 25;

% From your high-stim plot, there appear to be 7 pulses.
% If detection looks wrong, try setting this to [].
pulseParams.expectedNumPulses = 7;

pulseParams.maxNumPulses = 20;
pulseParams.thresholdMAD = 6;
pulseParams.localBaselineMS = 20;
pulseParams.templateStatistic = 'median';
pulseParams.taperEdgeMS = 2;

pulseParams.doBaselineCorrection = false;

pulseOut = PulsewiseTemplate(dataHighStim,timeVals,pulseParams);

fprintf('\nPulsewiseTemplate detected %d pulses.\n',length(pulseOut.pulseTimes));
fprintf('Detected pulse times, ms:\n');
disp(pulseOut.pulseTimes * 1000);

%% FFT summaries
fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);
fftERPClean = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);
fftAlignedClean = compute_fft_summary(erpAlignedOut.cleanedData,timeVals,fftIdx);
fftPulseClean = compute_fft_summary(pulseOut.cleanedData,timeVals,fftIdx);

%% Figure 1: Pulse detection and pulse-wise template
figure;

subplot(4,1,1);
plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.2);
hold on;
plot(pulseOut.pulseTimes,mean(dataHighStim(:,pulseOut.pulseSamples),1),'ko','MarkerFaceColor','k');
xlim([-0.1 0.5]);
title('High-stim mean with detected pulse locations');
xlabel('Time (s)');
ylabel('Mean LFP');

subplot(4,1,2);
plot(timeVals,pulseOut.detectionScore,'LineWidth',1.2);
hold on;
yline(pulseOut.detectionThreshold,'r--','LineWidth',1.2);
plot(pulseOut.pulseTimes,pulseOut.detectionScore(pulseOut.pulseSamples),'ko','MarkerFaceColor','k');
xlim([-0.1 0.5]);
title('Pulse detection score');
xlabel('Time (s)');
ylabel('Detection score');
legend('Detection score','Threshold','Detected pulses','Location','best');

subplot(4,1,3);
plot(timeVals,pulseOut.fullTemplate,'k','LineWidth',1.5);
xlim([-0.1 0.5]);
title('Full pulse-wise artifact template');
xlabel('Time (s)');
ylabel('Template');

subplot(4,1,4);
plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.1);
hold on;
plot(timeVals,pulseOut.fullTemplate,'k','LineWidth',1.3);
plot(timeVals,mean(pulseOut.cleanedData,1),'b','LineWidth',1.1);
xlim([-0.1 0.5]);
title('Mean signal before, pulse-wise template, and after correction');
xlabel('Time (s)');
ylabel('LFP');
legend('High-stim raw mean','Pulse-wise template','Cleaned mean','Location','best');

%% Figure 2: Before and after trial-level diagnostic
numTrialsToShow = 20;
trialSubset = 1:min(numTrialsToShow,size(dataHighStim,1));

figure;

subplot(3,1,1);
plot(timeVals,dataHighStim(trialSubset,:)');
xlim([-0.2 0.8]);
title('High-stim trials before PulsewiseTemplate correction');
xlabel('Time (s)');
ylabel('LFP');

subplot(3,1,2);
plot(timeVals,pulseOut.fullTemplate,'k','LineWidth',1.5);
xlim([-0.2 0.8]);
title('Pulse-wise template being subtracted');
xlabel('Time (s)');
ylabel('Template');

subplot(3,1,3);
plot(timeVals,pulseOut.cleanedData(trialSubset,:)');
xlim([-0.2 0.8]);
title('High-stim trials after PulsewiseTemplate correction');
xlabel('Time (s)');
ylabel('Cleaned LFP');

%% Figure 3: Individual pulse templates
figure;

imagesc(pulseOut.pulseOffsetTimes*1000,1:size(pulseOut.pulseTemplates,1),pulseOut.pulseTemplates);
colorbar;
title('Individual pulse templates');
xlabel('Time relative to detected pulse (ms)');
ylabel('Pulse number');

%% Figure 4: Mean time-series comparison
figure;

subplot(5,1,1);
plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.2);
xlim([-0.2 0.8]);
title('No-stim reference mean');
xlabel('Time (s)');
ylabel('LFP');

subplot(5,1,2);
plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.2);
xlim([-0.2 0.8]);
title('High-stim raw mean');
xlabel('Time (s)');
ylabel('LFP');

subplot(5,1,3);
plot(timeVals,mean(erpOut.cleanedData,1),'b','LineWidth',1.2);
xlim([-0.2 0.8]);
title('High-stim mean after ERPSubtraction');
xlabel('Time (s)');
ylabel('LFP');

subplot(5,1,4);
plot(timeVals,mean(erpAlignedOut.cleanedData,1),'g','LineWidth',1.2);
xlim([-0.2 0.8]);
title('High-stim mean after ERPAligned');
xlabel('Time (s)');
ylabel('LFP');

subplot(5,1,5);
plot(timeVals,mean(pulseOut.cleanedData,1),'m','LineWidth',1.2);
xlim([-0.2 0.8]);
title('High-stim mean after PulsewiseTemplate');
xlabel('Time (s)');
ylabel('LFP');

%% Figure 5: FFT comparison
figure;

plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.3);
hold on;
plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'r','LineWidth',1.2);
plot(fftERPClean.freqAxis,fftERPClean.logMeanMagnitude,'b','LineWidth',1.2);
plot(fftAlignedClean.freqAxis,fftAlignedClean.logMeanMagnitude,'g','LineWidth',1.3);
plot(fftPulseClean.freqAxis,fftPulseClean.logMeanMagnitude,'m','LineWidth',1.3);

xlim([0 200]);
title('FFT comparison: PulsewiseTemplate vs previous methods');
xlabel('Frequency (Hz)');
ylabel('log10 mean FFT magnitude');

legend( ...
    'No-stim reference', ...
    'High-stim raw', ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'PulsewiseTemplate', ...
    'Location','best');

%% Quantitative FFT metric
freqRangeForMetric = [0 200];

freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
           fftNoStim.freqAxis <= freqRangeForMetric(2);

rawError = norm(fftHighRaw.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

erpError = norm(fftERPClean.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

alignedError = norm(fftAlignedClean.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

pulseError = norm(fftPulseClean.logMeanMagnitude(freqMask) - ...
                  fftNoStim.logMeanMagnitude(freqMask));

erpImprovement = 100 * (rawError - erpError) / rawError;
alignedImprovement = 100 * (rawError - alignedError) / rawError;
pulseImprovement = 100 * (rawError - pulseError) / rawError;

fprintf('\nPulsewiseTemplate comparison results:\n');
fprintf('Raw FFT error vs no-stim          = %.4f\n',rawError);
fprintf('ERPSubtraction error              = %.4f\n',erpError);
fprintf('ERPAligned error                  = %.4f\n',alignedError);
fprintf('PulsewiseTemplate error           = %.4f\n',pulseError);

fprintf('\nImprovement relative to raw high-stim:\n');
fprintf('ERPSubtraction improvement        = %.2f %%\n',erpImprovement);
fprintf('ERPAligned improvement            = %.2f %%\n',alignedImprovement);
fprintf('PulsewiseTemplate improvement     = %.2f %%\n',pulseImprovement);