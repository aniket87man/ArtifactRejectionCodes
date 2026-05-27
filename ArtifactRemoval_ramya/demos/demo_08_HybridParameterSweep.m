clear; close all; clc;

%% demo_08_HybridParameterSweep
% Parameter sweep for the hybrid method:
% ERPAligned first, then PulsewiseTemplate on the ERPAligned residual.

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

%% Baseline methods
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);

alignedParams = struct();
alignedParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
alignedParams.alignWindow = [-0.01 0.03];
alignedParams.maxShiftMS = 10;
alignedParams.doBaselineCorrection = false;

erpAlignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);

%% FFT summaries for reference and baselines
fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);
fftERPClean = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);
fftAlignedClean = compute_fft_summary(erpAlignedOut.cleanedData,timeVals,fftIdx);

freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
           fftNoStim.freqAxis <= freqRangeForMetric(2);

rawError = norm(fftHighRaw.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

erpError = norm(fftERPClean.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

alignedError = norm(fftAlignedClean.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

erpImprovement = 100 * (rawError - erpError) / rawError;
alignedImprovement = 100 * (rawError - alignedError) / rawError;

fprintf('\nBaseline results:\n');
fprintf('Raw FFT error vs no-stim     = %.4f\n',rawError);
fprintf('ERPSubtraction error         = %.4f, improvement = %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned error             = %.4f, improvement = %.2f %%\n',alignedError,alignedImprovement);

%% Parameter sweep values
pulseWindowList = [
    -2  10
    -2  20
    -5  20
    -5  35
    -10 40
];

thresholdMADList = [2 3 4 6];

templateStatisticList = {'median','mean'};

%% Sweep
resultRows = [];

rowCounter = 0;

fprintf('\nRunning hybrid parameter sweep...\n');

for iWin = 1:size(pulseWindowList,1)

    for iThr = 1:length(thresholdMADList)

        for iStat = 1:length(templateStatisticList)

            pulseParams = struct();

            pulseParams.pulseSearchWindow = [-0.02 0.32];
            pulseParams.pulseWindowMS = pulseWindowList(iWin,:);
            pulseParams.minPulseDistanceMS = 25;
            pulseParams.expectedNumPulses = 7;
            pulseParams.maxNumPulses = 20;
            pulseParams.thresholdMAD = thresholdMADList(iThr);
            pulseParams.localBaselineMS = 20;
            pulseParams.templateStatistic = templateStatisticList{iStat};
            pulseParams.taperEdgeMS = 2;
            pulseParams.doBaselineCorrection = false;

            try
                hybridParams = struct();
                hybridParams.alignedParams = alignedParams;
                hybridParams.pulseParams = pulseParams;

                hybridOut = ERPAlignedPulsewise(dataHighStim,timeVals,hybridParams);

                fftHybridClean = compute_fft_summary(hybridOut.cleanedData,timeVals,fftIdx);

                hybridError = norm(fftHybridClean.logMeanMagnitude(freqMask) - ...
                                   fftNoStim.logMeanMagnitude(freqMask));

                hybridImprovement = 100 * (rawError - hybridError) / rawError;

                numPulsesDetected = length(hybridOut.pulseOut.pulseTimes);
                clippedTrials = sum(hybridOut.erpAlignedOut.wasClipped);

                rowCounter = rowCounter + 1;

                resultRows(rowCounter).pulseWindowStartMS = pulseWindowList(iWin,1);
                resultRows(rowCounter).pulseWindowEndMS = pulseWindowList(iWin,2);
                resultRows(rowCounter).thresholdMAD = thresholdMADList(iThr);
                resultRows(rowCounter).templateStatistic = string(templateStatisticList{iStat});
                resultRows(rowCounter).numPulsesDetected = numPulsesDetected;
                resultRows(rowCounter).hybridError = hybridError;
                resultRows(rowCounter).hybridImprovement = hybridImprovement;
                resultRows(rowCounter).clippedTrials = clippedTrials;

                fprintf('Done: window [%d %d] ms, threshold %.1f, %s -> error %.4f, improvement %.2f %%\n', ...
                    pulseWindowList(iWin,1),pulseWindowList(iWin,2), ...
                    thresholdMADList(iThr),templateStatisticList{iStat}, ...
                    hybridError,hybridImprovement);

            catch ME

                rowCounter = rowCounter + 1;

                resultRows(rowCounter).pulseWindowStartMS = pulseWindowList(iWin,1);
                resultRows(rowCounter).pulseWindowEndMS = pulseWindowList(iWin,2);
                resultRows(rowCounter).thresholdMAD = thresholdMADList(iThr);
                resultRows(rowCounter).templateStatistic = string(templateStatisticList{iStat});
                resultRows(rowCounter).numPulsesDetected = NaN;
                resultRows(rowCounter).hybridError = NaN;
                resultRows(rowCounter).hybridImprovement = NaN;
                resultRows(rowCounter).clippedTrials = NaN;

                fprintf('FAILED: window [%d %d] ms, threshold %.1f, %s\n', ...
                    pulseWindowList(iWin,1),pulseWindowList(iWin,2), ...
                    thresholdMADList(iThr),templateStatisticList{iStat});

                fprintf('Reason: %s\n',ME.message);

            end
        end
    end
end

%% Convert to table and sort
resultsTable = struct2table(resultRows);

validRows = ~isnan(resultsTable.hybridError);
validResults = resultsTable(validRows,:);

validResults = sortrows(validResults,'hybridError','ascend');

fprintf('\nTop hybrid parameter settings:\n');
disp(validResults(1:min(10,height(validResults)),:));

%% Best setting
bestRow = validResults(1,:);

fprintf('\nBest hybrid setting:\n');
disp(bestRow);

%% Re-run best hybrid for plotting
bestPulseParams = struct();
bestPulseParams.pulseSearchWindow = [-0.02 0.32];
bestPulseParams.pulseWindowMS = [bestRow.pulseWindowStartMS bestRow.pulseWindowEndMS];
bestPulseParams.minPulseDistanceMS = 25;
bestPulseParams.expectedNumPulses = 7;
bestPulseParams.maxNumPulses = 20;
bestPulseParams.thresholdMAD = bestRow.thresholdMAD;
bestPulseParams.localBaselineMS = 20;
bestPulseParams.templateStatistic = char(bestRow.templateStatistic);
bestPulseParams.taperEdgeMS = 2;
bestPulseParams.doBaselineCorrection = false;

bestHybridParams = struct();
bestHybridParams.alignedParams = alignedParams;
bestHybridParams.pulseParams = bestPulseParams;

bestHybridOut = ERPAlignedPulsewise(dataHighStim,timeVals,bestHybridParams);

fftBestHybridClean = compute_fft_summary(bestHybridOut.cleanedData,timeVals,fftIdx);

%% Figure 1: parameter sweep results
figure;
scatter(validResults.hybridError,validResults.hybridImprovement,60,'filled');
xlabel('FFT error vs no-stim');
ylabel('Improvement relative to raw high-stim (%)');
title('Hybrid parameter sweep results');
grid on;

%% Figure 2: top 10 settings
topN = min(10,height(validResults));

figure;
bar(validResults.hybridImprovement(1:topN));
xticks(1:topN);
xlabel('Top parameter setting rank');
ylabel('Improvement (%)');
title('Top hybrid parameter settings');

%% Figure 3: FFT comparison using best hybrid setting
figure;

plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.3);
hold on;
plot(fftHighRaw.freqAxis,fftHighRaw.logMeanMagnitude,'r','LineWidth',1.2);
plot(fftERPClean.freqAxis,fftERPClean.logMeanMagnitude,'b','LineWidth',1.2);
plot(fftAlignedClean.freqAxis,fftAlignedClean.logMeanMagnitude,'g','LineWidth',1.3);
plot(fftBestHybridClean.freqAxis,fftBestHybridClean.logMeanMagnitude,'m','LineWidth',1.3);

xlim([0 200]);
title('FFT comparison using best hybrid parameter setting');
xlabel('Frequency (Hz)');
ylabel('log10 mean FFT magnitude');

legend( ...
    'No-stim reference', ...
    'High-stim raw', ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'Best ERPAlignedPulsewise', ...
    'Location','best');

%% Figure 4: best hybrid time-domain diagnostic
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
plot(timeVals,bestHybridOut.erpAlignedOut.cleanedData(trialSubset,:)');
xlim([-0.2 0.8]);
title('After ERPAligned stage');
xlabel('Time (s)');
ylabel('LFP');

subplot(4,1,3);
plot(timeVals,bestHybridOut.pulseOut.fullTemplate,'k','LineWidth',1.5);
xlim([-0.2 0.8]);
title('Best residual pulse-wise template');
xlabel('Time (s)');
ylabel('Template');

subplot(4,1,4);
plot(timeVals,bestHybridOut.cleanedData(trialSubset,:)');
xlim([-0.2 0.8]);
title('After best ERPAlignedPulsewise hybrid');
xlabel('Time (s)');
ylabel('Cleaned LFP');

%% Print final summary
bestError = bestRow.hybridError;
bestImprovement = bestRow.hybridImprovement;

fprintf('\nFinal best result from sweep:\n');
fprintf('Raw FFT error              = %.4f\n',rawError);
fprintf('ERPSubtraction error       = %.4f, improvement = %.2f %%\n',erpError,erpImprovement);
fprintf('ERPAligned error           = %.4f, improvement = %.2f %%\n',alignedError,alignedImprovement);
fprintf('Best hybrid error          = %.4f, improvement = %.2f %%\n',bestError,bestImprovement);

fprintf('\nBest parameters:\n');
fprintf('pulseWindowMS              = [%g %g]\n',bestRow.pulseWindowStartMS,bestRow.pulseWindowEndMS);
fprintf('thresholdMAD               = %.2f\n',bestRow.thresholdMAD);
fprintf('templateStatistic          = %s\n',bestRow.templateStatistic);
fprintf('numPulsesDetected          = %d\n',bestRow.numPulsesDetected);