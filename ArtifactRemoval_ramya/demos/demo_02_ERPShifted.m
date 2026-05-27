%% demo_02_ERPShifted
% ERPShifted comparison with ERPSubtraction baseline
% demo for elec1 and elec7
%
% Methods compared:
%   1. No-stimulation reference
%   2. Uncleaned high-stimulation signal
%   3. ERPSubtraction
%   4. ERPShifted
%
% Important:
% ERPShifted.m should estimate shifts using the full ERP template,
% but subtract the shifted template only inside the artifact window.

clear; close all; clc;

%% Paths
folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';
subjectName = 'dona';
gridType = 'Microelectrode';
expDate = '290825';
protocolName = 'GRF_001';

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

figFolder = fullfile(pwd,'figures');
if ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% Load common files
I = load(lfpInfoFile);
P = load(paramFile);

timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

%% Conditions
noStimTrials = parameterCombinations{1,1,1,5,5,4};
highStimTrials = parameterCombinations{7,1,1,5,5,4};

%% Windows
artifactWindow = [0 0.4];

% For clean time-domain plots, stop at 0.4 s because subtraction is only
% applied inside this window.
plotWindow = [-0.05 0.4];

% FFT window.
% This gives 800 samples: 0.0005 to 0.4000 s
pos = find(timeVals > artifactWindow(1) & timeVals <= artifactWindow(2));

Fs = 1/median(diff(timeVals));
N = length(pos);
freqVals = (0:N-1) * (Fs/N);

freqRangeForMetric = [0 200];
freqMask = freqVals >= freqRangeForMetric(1) & freqVals <= freqRangeForMetric(2);

fprintf('\n============================================================\n');
fprintf('demo_02_ERPShifted\n');
fprintf('============================================================\n');
fprintf('Artifact/subtraction window: %.2f to %.2f s\n',artifactWindow(1),artifactWindow(2));
fprintf('Actual FFT window: %.4f to %.4f s\n',timeVals(pos(1)),timeVals(pos(end)));
fprintf('FFT samples: %d\n',N);
fprintf('Frequency resolution: %.2f Hz\n',Fs/N);
fprintf('Frequency metric range: %.1f to %.1f Hz\n',freqRangeForMetric(1),freqRangeForMetric(2));

%% Electrodes to run
% elec1: stimulation-electrode illustration
% elec7: example analyzed electrode
electrodesToRun = {'elec1','elec7'};

summaryRows = {};

for iElec = 1:length(electrodesToRun)

    electrodeName = electrodesToRun{iElec};

    fprintf('\n============================================================\n');
    fprintf('Running ERPShifted comparison for %s\n',electrodeName);
    fprintf('============================================================\n');

    %% Load electrode data
    lfpFile = fullfile(lfpFolder,[electrodeName '.mat']);
    D = load(lfpFile);

    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    fprintf('No-stim trials   : %d\n',size(dataNoStim,1));
    fprintf('High-stim trials : %d\n',size(dataHighStim,1));

    %% Method 1: ERPSubtraction baseline
    erpParams = struct();
    erpParams.subtractWindow = artifactWindow;
    erpParams.doBaselineCorrection = false;

    erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);

    %% Method 2: ERPShifted
    shiftParams = struct();
    shiftParams.subtractWindow = artifactWindow;
    shiftParams.alignWindow = [-0.01 0.03];
    shiftParams.maxShiftMS = 10;
    shiftParams.doBaselineCorrection = false;

    erpShiftOut = ERPShifted(dataHighStim,timeVals,shiftParams);

    %% Sanity checks
    subtractIdx = timeVals > artifactWindow(1) & timeVals <= artifactWindow(2);
    outsideIdx = ~subtractIdx;

    erpOutsideChange = max(abs(erpOut.cleanedData(:,outsideIdx) - ...
                               dataHighStim(:,outsideIdx)),[],'all');

    shiftedOutsideChange = max(abs(erpShiftOut.cleanedData(:,outsideIdx) - ...
                                   dataHighStim(:,outsideIdx)),[],'all');

    fprintf('\nSanity checks:\n');
    fprintf('ERPSubtraction max outside-window change = %.12f\n',erpOutsideChange);
    fprintf('ERPShifted max outside-window change     = %.12f\n',shiftedOutsideChange);

    %% FFT calculation
    fftNoStim = log10(mean(abs(fft(dataNoStim(:,pos)'))'));
    fftHighRaw = log10(mean(abs(fft(dataHighStim(:,pos)'))'));
    fftERPClean = log10(mean(abs(fft(erpOut.cleanedData(:,pos)'))'));
    fftShiftClean = log10(mean(abs(fft(erpShiftOut.cleanedData(:,pos)'))'));

    %% Quantitative FFT error
    rawError = norm(fftHighRaw(freqMask) - fftNoStim(freqMask));
    erpError = norm(fftERPClean(freqMask) - fftNoStim(freqMask));
    shiftError = norm(fftShiftClean(freqMask) - fftNoStim(freqMask));

    erpImprovement = 100 * (rawError - erpError) / rawError;
    shiftImprovement = 100 * (rawError - shiftError) / rawError;

    fprintf('\nFFT comparison results for %s:\n',electrodeName);
    fprintf('Raw FFT error vs no-stim         = %.4f\n',rawError);
    fprintf('ERPSubtraction error             = %.4f\n',erpError);
    fprintf('ERPShifted error                 = %.4f\n',shiftError);
    fprintf('ERPSubtraction improvement       = %.2f %%\n',erpImprovement);
    fprintf('ERPShifted improvement           = %.2f %%\n',shiftImprovement);

    %% Shift summary
    shiftMS = erpShiftOut.shiftTimes * 1000;

    minShiftMS = min(shiftMS);
    maxShiftMS = max(shiftMS);
    meanShiftMS = mean(shiftMS);
    stdShiftMS = std(shiftMS);

    fprintf('\nShift summary for %s:\n',electrodeName);
    fprintf('Min shift  = %d samples, %.4f ms\n', ...
        min(erpShiftOut.shiftSamples),minShiftMS);
    fprintf('Max shift  = %d samples, %.4f ms\n', ...
        max(erpShiftOut.shiftSamples),maxShiftMS);
    fprintf('Mean shift = %.2f samples, %.4f ms\n', ...
        mean(erpShiftOut.shiftSamples),meanShiftMS);
    fprintf('Std shift  = %.2f samples, %.4f ms\n', ...
        std(erpShiftOut.shiftSamples),stdShiftMS);

    summaryRows(end+1,:) = {electrodeName, ...
        rawError,erpError,shiftError, ...
        erpImprovement,shiftImprovement, ...
        minShiftMS,maxShiftMS,meanShiftMS,stdShiftMS, ...
        erpOutsideChange,shiftedOutsideChange};

    %% Figure 1: Time-domain diagnostic
    numTrialsToShow = 30;
    trialSubset = 1:min(numTrialsToShow,size(dataHighStim,1));

    figure('Color','w','Position',[100 100 1100 850]);

    subplot(4,1,1);
    plot(timeVals,dataHighStim(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(dataHighStim,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['Uncleaned high-stimulation trials: ' electrodeName]);
    box off;

    subplot(4,1,2);
    plot(timeVals,erpOut.template,'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('ERP template');
    box off;

    subplot(4,1,3);
    plot(timeVals,erpShiftOut.shiftedTemplates(trialSubset,:)','Color',[0.70 0.70 0.70]);
    hold on;
    plot(timeVals,mean(erpShiftOut.shiftedTemplates,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('Trial-wise shifted ERP templates');
    box off;

    subplot(4,1,4);
    plot(timeVals,erpShiftOut.cleanedData(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(erpShiftOut.cleanedData,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['ERPShifted residual trials: ' electrodeName]);
    box off;

    sgtitle(['ERPShifted template subtraction: ' electrodeName]);

    outFileTime = fullfile(figFolder,['demo02_shifted_time_' electrodeName '.png']);
    exportgraphics(gcf,outFileTime,'Resolution',300);

    %% Figure 2: FFT comparison
    figure('Color','w','Position',[100 100 900 550]);

    plot(freqVals,fftNoStim,'k','LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'r','LineWidth',1.7);
    plot(freqVals,fftERPClean,'b','LineWidth',1.7);
    plot(freqVals,fftShiftClean,'m','LineWidth',1.7);

    xlim([0 200]);

    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');
    title(['FFT comparison: ERPSubtraction vs ERPShifted, ' electrodeName]);

    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERPSubtraction', ...
            'ERPShifted'}, ...
            'Location','best');

    box off;

    outFileFFT = fullfile(figFolder,['demo02_shifted_fft_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileFFT,'ContentType','vector');

    %% Figure 3: Shift distribution
    figure('Color','w','Position',[100 100 750 450]);

    histogram(shiftMS,'BinWidth',0.5);
    hold on;
    xline(0,'k--','LineWidth',1.2);
    xline(meanShiftMS,'r--','LineWidth',1.2);

    title(['Distribution of estimated ERP template shifts: ' electrodeName]);
    xlabel('Estimated shift (ms)');
    ylabel('Number of trials');
    legend({'Shift distribution','Zero shift','Mean shift'},'Location','best');
    box off;

    outFileShift = fullfile(figFolder,['demo02_shifted_shift_hist_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileShift,'ContentType','vector');

    fprintf('\nSaved:\n');
    fprintf('%s\n',outFileTime);
    fprintf('%s\n',outFileFFT);
    fprintf('%s\n',outFileShift);

end

%% Save summary table
summaryTable = cell2table(summaryRows, ...
    'VariableNames',{'Electrode', ...
    'RawFFTError','ERPFFTError','ShiftedFFTError', ...
    'ERPImprovementPercent','ShiftedImprovementPercent', ...
    'MinShiftMS','MaxShiftMS','MeanShiftMS','StdShiftMS', ...
    'ERPOutsideWindowChange','ShiftedOutsideWindowChange'});

disp(summaryTable);

summaryFile = fullfile(figFolder,'demo02_shifted_summary.csv');
writetable(summaryTable,summaryFile);

fprintf('\nSaved ERPShifted summary:\n%s\n',summaryFile);
fprintf('\ndemo_02_ERPShifted complete.\n');