clear; 
% close all;

subjectName = 'dona';

expDate = '290825'; protocolName = 'GRF_001'; folderSourceString = '/Users/aniketmandal/Documents/MATLAB/SupiLabProgramsDatas'; 
gridType = 'Microelectrode'; stimElectrode = 1; numElec = 36;

lfpData = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(stimElectrode) '.mat']));
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', 'lfpInfo.mat'));
timeVals = tmp.timeVals;

tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'extractedData', 'parameterCombinations.mat'));
p = tmp.parameterCombinations;

Fs = round(1/(timeVals(2) - timeVals(1)));
stRange = [0 0.4];
pos = intersect(find(timeVals>stRange(1)), find(timeVals<stRange(2)));      % stimPos

df = 1/diff(stRange);
freqVals = 0: df : (Fs - df);

a = 7;                      % microStim amp.
dataNoStim = lfpData.analogData(p{1,1,1,5,5,4},pos);
dataStim = lfpData.analogData(p{a,1,1,5,5,4},pos);

figure;
subplot(121)
plot(timeVals(pos),mean(dataNoStim,1),'g'); hold on;
plot(timeVals(pos),mean(dataStim,1),'r');


subplot(122)
plot(freqVals,log10(mean(abs(fft(dataNoStim')'))),'g'); hold on;
plot(freqVals,log10(mean(abs(fft(dataStim')'))),'k');

% Subtract ERP
dataNoStimMeanSubracted = dataNoStim - repmat(mean(dataNoStim,1),size(dataNoStim,1),1);
dataStimMeanSubracted = dataStim - repmat(mean(dataStim,1),size(dataStim,1),1);

figure;
plot(freqVals,log10(mean(abs(fft(dataNoStimMeanSubracted')'))),'c'); hold on;
plot(freqVals,log10(mean(abs(fft(dataStimMeanSubracted')'))),'r');
xlim([0 500])

% % % % % % 
% test stimStartPosn aligning
fullDataStim = lfpData.analogData(p{a,1,1,5,5,4},:);
[dataStimStartPosns, alignedData] = alignStimStartPosns(fullDataStim, timeVals);

alignedDataStim = alignedData(:, pos);     % only looking at stimPos

% plot and check how well it does againt original, non-aligned Stim data
figure;
subplot(121)
plot(timeVals(pos),mean(dataNoStim,1),'g'); hold on;
plot(timeVals(pos),mean(dataStim,1),'r'); plot(timeVals(pos), mean(alignedDataStim), 'k');
legend('No microStim.', 'with microStim.', 'microStim data aligned to stimStartPosn.');

subplot(122)
plot(freqVals,log10(mean(abs(fft(dataNoStim')'))),'g'); hold on;
plot(freqVals,log10(mean(abs(fft(dataStim')'))),'r');
plot(freqVals,log10(mean(abs(fft(alignedDataStim')'))),'k');     % looks higher cause all the peaks more aligned,
legend('No microStim.', 'with microStim.', 'microStim data aligned to stimStartPosn.');   % so more averaging in the non-aligned case
                                                                 
% subtract ERP
alignedDataStimMeanSubracted = alignedDataStim - repmat(mean(alignedDataStim,1),size(alignedDataStim,1),1);

figure;      % plot against original ERP subtracted PSD
plot(freqVals,log10(mean(abs(fft(dataNoStimMeanSubracted')'))),'c'); hold on;
plot(freqVals,log10(mean(abs(fft(dataStimMeanSubracted')'))),'r');
plot(freqVals,log10(mean(abs(fft(alignedDataStimMeanSubracted')'))),'k');
xlim([0 500])
title('PSD Comparison'); legend('No microStim.', 'with microStim.', 'microStim data aligned to stimStartPosn.');

figure;
plot(timeVals, fullDataStim - mean(fullDataStim, 1)); title('non-aligned Stimmed Data, Mean-subtracted')

figure;
plot(timeVals(pos), alignedDataStimMeanSubracted); title('aligned Stimmed Data, Mean-subtracted')

figure;
plot(timeVals, alignedData - mean(alignedData, 1)); title('aligned Stimmed Data, Mean-subtracted')


% apply optimal jitter, scale to ERP before subtracting
alignedDataMean = mean(alignedData, 1);
[processedAlignedData, ~] = applyOptimalJitterScale(alignedData, alignedDataMean, []);

% plot data after subtracting optimally jittered and scaled ERP from the alignedData
figure;
plot(timeVals, processedAlignedData);
title('aligned Stimmed Data, Mean-subtracted, optimal Jitter & Scale applied')
% figure(10); hold on; plot(timeVals, alignedDataMean, 'k')
% xlim([-0.01 0.4])

% processedAlignedDataStim = applyOptimalJitterScale(alignedDataStim, mean(alignedDataStim, 1));   % alternate way of getting processedAlignedDataStim

% figure;
% plot(timeVals(pos), processedAlignedDataStim);
% title('aligned Stimmed Data, Mean-subtracted, optimal Jitter & Scale applied between 0 and 0.4s')

processedAlignedDataStim = processedAlignedData(:, pos);

figure;      % plot against original ERP subtracted PSD
plot(freqVals,log10(mean(abs(fft(dataNoStimMeanSubracted')'))),'c'); hold on;
plot(freqVals,log10(mean(abs(fft(dataStimMeanSubracted')'))),'r');
plot(freqVals,log10(mean(abs(fft(alignedDataStimMeanSubracted')'))),'k');
plot(freqVals,log10(mean(abs(fft(processedAlignedDataStim')'))),'b');
xlim([0 500])
title('PSD Comparison'); legend('No microStim.', 'with microStim.', 'microStim data aligned to stimStartPosn.', 'alignedData, optimal Jitter & Scale applied ERP subtracted');

% Something might be wrong in my code


% apply optimal jitter, scale to ERP, based on only the spikes part, before subtracting
fullDataStim;
alignedData;
lastSpikePosn = find(alignedDataMean == max(alignedDataMean(timeVals>0.3 & timeVals<0.38)));

jitterWindow = timeVals>-0.5 & timeVals<(timeVals(lastSpikePosn)+0.01);
[alignedDataSpikeJittered, jitteredScaledErps] = applyOptimalJitterScale(alignedData, alignedDataMean, jitterWindow);

% plot data after subtracting optimally jittered and scaled ERP from the alignedData
figure;
plot(timeVals, alignedDataSpikeJittered);
title('aligned Stimmed Data, Mean-subtracted, optimal Jitter & Scale applied according to only the artifact spikes')

    % try on the unaligned data
fullDataJittered = applyOptimalJitterScale(fullDataStim, mean(fullDataStim, 1), []);
figure;
plot(timeVals, fullDataJittered); title('unaligned Data Jittered & Scaled')


% other electrodes
elec2_1 = 7; elec2_2 = 8;   % electrodes in the 2nd distance grp
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(elec2_1) '.mat']));
fullDataElec2_1 = tmp.analogData(p{a, 1,1,5,5,4}, :);
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(elec2_2) '.mat']));
fullDataElec2_2 = tmp.analogData(p{a, 1,1,5,5,4}, :);

% try on the unaligned data of elec7 and elec8
fullDataJitteredElec2_1 = applyOptimalJitterScale(fullDataElec2_1, mean(fullDataElec2_1, 1), jitterWindow);
figure;
plot(timeVals, fullDataJitteredElec2_1); title(sprintf('elec%d Data Jittered & Scaled based on spike Window', elec2_1))


elec5_1 = 22;     % electrodes in the 5th distance group
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(elec5_1) '.mat']));
fullDataElec5_1 = tmp.analogData(p{a, 1,1,5,5,4}, :);

% try on the unaligned data of elec22
fullDataJitteredElec5_1 = applyOptimalJitterScale(fullDataElec5_1, mean(fullDataElec5_1, 1), []);
figure;
plot(timeVals, fullDataJitteredElec5_1); title(sprintf('elec%d Data Jittered & Scaled', elec5_1))

[optimShift_lagData_1_7, maxCorr_1_7] = getCorrelations(optimalShifts_unAligned', optimalShifts_elec7');
[optimShift_lagData_1_22, maxCorr_1_22] = getCorrelations(optimalShifts_unAligned', optimalShifts_elec22');
[optimShift_lagData_7_22, maxCorr_7_22] = getCorrelations(optimalShifts_elec7', optimalShifts_elec22');

% % % % % on other microStimn. amplitudes % % % % %
a=3;
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(stimElectrode) '.mat']));
fullDataElec1 = tmp.analogData(p{a,1,1,5,5,4}, :);

fullDataJitteredElec1 = applyOptimalJitterScale(fullDataElec1, mean(fullDataElec1, 1), jitterWindow);
figure;
plot(timeVals, fullDataJitteredElec1); title(sprintf('elec%d Data Jittered & Scaled', stimElectrode))

elec2_1 = 7; elec2_2 = 8;   % electrodes in the 2nd distance grp
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(elec2_1) '.mat']));
fullDataElec2_1 = tmp.analogData(p{a,1,1,5,5,4}, :);
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(elec2_2) '.mat']));
fullDataElec2_2 = tmp.analogData(p{a,1,1,5,5,4}, :);

% try on the unaligned data of elec7 and elec8
fullDataJitteredElec2_1 = applyOptimalJitterScale(fullDataElec2_1, mean(fullDataElec2_1, 1), []);
figure;
plot(timeVals, fullDataJitteredElec2_1); title(sprintf('elec%d Data Jittered & Scaled', elec2_1))


elec5_1 = 22;     % electrodes in the 5th distance group
tmp = load(fullfile(folderSourceString, 'data', subjectName, gridType, expDate, protocolName, 'segmentedData', 'LFP', ['elec' num2str(elec5_1) '.mat']));
fullDataElec5_1 = tmp.analogData(p{a, 1,1,5,5,4}, :);

% try on the unaligned data of elec22
fullDataJitteredElec5_1 = applyOptimalJitterScale(fullDataElec5_1, mean(fullDataElec5_1, 1), []);
figure;
plot(timeVals, fullDataJitteredElec5_1); title(sprintf('elec%d Data Jittered & Scaled based on spike Window', elec5_1))

[Shift_lagData_1_7, Corr_1_7] = getCorrelations(Shifts_elec1', Shifts_elec7');
[Shift_lagData_1_22, Corr_1_22] = getCorrelations(Shifts_elec1', Shifts_elec22');
[Shift_lagData_7_22, Corr_7_22] = getCorrelations(Shifts_elec7', Shifts_elec22');
