%% demo_05_ERPAligned
% Compare ERP-based variants including ERPAligned:
%   1. No-stimulation reference
%   2. Uncleaned high-stimulation signal
%   3. ERPSubtraction
%   4. ERPShifted
%   5. ERPScaled
%   6. ERPShiftedScaled
%   7. ERPAligned
%e 
% Runs elec1 and elec7.
%
% Final window convention:
%   [0 0.4] means timeVals > 0 and timeVals <= 0.4
%   Actual FFT window: 0.0005 to 0.4000 s
%   FFT samples: 800

clear; close all; clc;

clear ERPSubtraction ERPShifted ERPScaled ERPShiftedScaled ERPAligned

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
plotWindow = [-0.05 0.4];

% FFT window: 0.0005 to 0.4000 s, 800 samples
pos = find(timeVals > artifactWindow(1) & timeVals <= artifactWindow(2));

Fs = 1/median(diff(timeVals));
N = length(pos);
freqVals = (0:N-1) * (Fs/N);

freqRangeForMetric = [0 200];
freqMask = freqVals >= freqRangeForMetric(1) & freqVals <= freqRangeForMetric(2);

fprintf('\n============================================================\n');
fprintf('demo_05_ERPAligned\n');
fprintf('============================================================\n');
fprintf('Artifact/subtraction window: %.2f to %.2f s\n',artifactWindow(1),artifactWindow(2));
fprintf('Actual FFT window: %.4f to %.4f s\n',timeVals(pos(1)),timeVals(pos(end)));
fprintf('FFT samples: %d\n',N);
fprintf('Frequency resolution: %.2f Hz\n',Fs/N);
fprintf('Frequency metric range: %.1f to %.1f Hz\n',freqRangeForMetric(1),freqRangeForMetric(2));

%% Electrodes to run
electrodesToRun = {'elec1','elec7'};

summaryRows = {};

for iElec = 1:length(electrodesToRun)

    electrodeName = electrodesToRun{iElec};

    fprintf('\n============================================================\n');
    fprintf('Running ERPAligned comparison for %s\n',electrodeName);
    fprintf('============================================================\n');

    %% Load electrode data
    lfpFile = fullfile(lfpFolder,[electrodeName '.mat']);
    D = load(lfpFile);

    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    fprintf('No-stim trials   : %d\n',size(dataNoStim,1));
    fprintf('High-stim trials : %d\n',size(dataHighStim,1));

    %% Method 1: ERPSubtraction
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

    %% Method 3: ERPScaled
    scaleParams = struct();
    scaleParams.subtractWindow = artifactWindow;
    scaleParams.scaleWindow = artifactWindow;
    scaleParams.scaleBounds = [0 2];
    scaleParams.doBaselineCorrection = false;

    erpScaleOut = ERPScaled(dataHighStim,timeVals,scaleParams);

    %% Method 4: ERPShiftedScaled
    shiftScaleParams = struct();
    shiftScaleParams.subtractWindow = artifactWindow;
    shiftScaleParams.fitWindow = [-0.01 0.03];
    shiftScaleParams.maxShiftMS = 10;
    shiftScaleParams.scaleBounds = [0 2];
    shiftScaleParams.doBaselineCorrection = false;

    erpShiftScaleOut = ERPShiftedScaled(dataHighStim,timeVals,shiftScaleParams);

    %% Method 5: ERPAligned
    alignedParams = struct();
    alignedParams.subtractWindow = artifactWindow;
    alignedParams.alignWindow = [-0.01 0.03];
    alignedParams.maxShiftMS = 10;
    alignedParams.doBaselineCorrection = false;

    erpAlignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);

    %% Sanity checks
    subtractIdx = timeVals > artifactWindow(1) & timeVals <= artifactWindow(2);
    outsideIdx = ~subtractIdx;

    erpOutsideChange = max(abs(erpOut.cleanedData(:,outsideIdx) - ...
                               dataHighStim(:,outsideIdx)),[],'all');

    shiftedOutsideChange = max(abs(erpShiftOut.cleanedData(:,outsideIdx) - ...
                                   dataHighStim(:,outsideIdx)),[],'all');

    scaledOutsideChange = max(abs(erpScaleOut.cleanedData(:,outsideIdx) - ...
                                  dataHighStim(:,outsideIdx)),[],'all');

    shiftScaleOutsideChange = max(abs(erpShiftScaleOut.cleanedData(:,outsideIdx) - ...
                                      dataHighStim(:,outsideIdx)),[],'all');

    alignedOutsideChange = max(abs(erpAlignedOut.cleanedData(:,outsideIdx) - ...
                                   dataHighStim(:,outsideIdx)),[],'all');

    fprintf('\nSanity checks:\n');
    fprintf('ERPSubtraction max outside-window change    = %.12f\n',erpOutsideChange);
    fprintf('ERPShifted max outside-window change        = %.12f\n',shiftedOutsideChange);
    fprintf('ERPScaled max outside-window change         = %.12f\n',scaledOutsideChange);
    fprintf('ERPShiftedScaled max outside-window change  = %.12f\n',shiftScaleOutsideChange);
    fprintf('ERPAligned max outside-window change        = %.12f\n',alignedOutsideChange);

    %% FFT calculation using trial-wise magniude averaging
    fftNoStim = log10(mean(abs(fft(dataNoStim(:,pos)'))'));
    fftHighRaw = log10(mean(abs(fft(dataHighStim(:,pos)'))'));
    fftERP = log10(mean(abs(fft(erpOut.cleanedData(:,pos)'))'));
    fftShift = log10(mean(abs(fft(erpShiftOut.cleanedData(:,pos)'))'));
    fftScale = log10(mean(abs(fft(erpScaleOut.cleanedData(:,pos)'))'));
    fftShiftScale = log10(mean(abs(fft(erpShiftScaleOut.cleanedData(:,pos)'))'));
    fftAligned = log10(mean(abs(fft(erpAlignedOut.cleanedData(:,pos)'))'));


    %% Figure: Comparison between ERPSubtraction and ERPAligned

    figure('Color','w','Position',[100 100 1100 750]);

    % Panel 1: Time-domain mean comparison
    subplot(2,1,1);

    plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.7);
    hold on;
    plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.5);
    plot(timeVals,mean(erpOut.cleanedData,1),'b','LineWidth',1.5);
    plot(timeVals,mean(erpAlignedOut.cleanedData,1),'g','LineWidth',1.7);

    xlim(plotWindow);

    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['Time-domain comparison: ERPSubtraction vs ERPAligned, ' electrodeName]);

    legend({'No-stimulation reference', ...
        'Uncleaned 64 \muA', ...
        'ERPSubtraction', ...
        'ERPAligned'}, ...
        'Location','best');

    box off;

    % Panel 2: FFT comparison
    subplot(2,1,2);

    plot(freqVals,fftNoStim,'k','LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'r','LineWidth',1.5);
    plot(freqVals,fftERP,'b','LineWidth',1.5);
    plot(freqVals,fftAligned,'g','LineWidth',1.7);

    xlim([0 200]);

    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');
    title(['Frequency-domain comparison: ERPSubtraction vs ERPAligned, ' electrodeName]);

    legend({'No-stimulation reference', ...
        'Uncleaned 64 \muA', ...
        'ERPSubtraction', ...
        'ERPAligned'}, ...
        'Location','best');

    box off;

    sgtitle(['ERPAligned comparison with ERPSubtraction baseline: ' electrodeName]);

    outFileSimple = fullfile(figFolder,['demo05_simple_erp_vs_aligned_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileSimple,'ContentType','vector');

    fprintf('%s\n',outFileSimple);


    %% Quantitative FFT errors
    rawError = norm(fftHighRaw(freqMask) - fftNoStim(freqMask));
    erpError = norm(fftERP(freqMask) - fftNoStim(freqMask));
    shiftError = norm(fftShift(freqMask) - fftNoStim(freqMask));
    scaleError = norm(fftScale(freqMask) - fftNoStim(freqMask));
    shiftScaleError = norm(fftShiftScale(freqMask) - fftNoStim(freqMask));
    alignedError = norm(fftAligned(freqMask) - fftNoStim(freqMask));

    erpImprovement = 100 * (rawError - erpError) / rawError;
    shiftImprovement = 100 * (rawError - shiftError) / rawError;
    scaleImprovement = 100 * (rawError - scaleError) / rawError;
    shiftScaleImprovement = 100 * (rawError - shiftScaleError) / rawError;
    alignedImprovement = 100 * (rawError - alignedError) / rawError;

    fprintf('\nFFT comparison results for %s:\n',electrodeName);
    fprintf('Raw FFT error vs no-stim         = %.4f\n',rawError);
    fprintf('ERPSubtraction error             = %.4f\n',erpError);
    fprintf('ERPShifted error                 = %.4f\n',shiftError);
    fprintf('ERPScaled error                  = %.4f\n',scaleError);
    fprintf('ERPShiftedScaled error           = %.4f\n',shiftScaleError);
    fprintf('ERPAligned error                 = %.4f\n',alignedError);

    fprintf('\nImprovement over uncleaned high-stim:\n');
    fprintf('ERPSubtraction improvement       = %.2f %%\n',erpImprovement);
    fprintf('ERPShifted improvement           = %.2f %%\n',shiftImprovement);
    fprintf('ERPScaled improvement            = %.2f %%\n',scaleImprovement);
    fprintf('ERPShiftedScaled improvement     = %.2f %%\n',shiftScaleImprovement);
    fprintf('ERPAligned improvement           = %.2f %%\n',alignedImprovement);

    %% ERPAligned shift summary
    alignedShiftMS = erpAlignedOut.shiftTimes * 1000;

    alignedMinShiftMS = min(alignedShiftMS);
    alignedMaxShiftMS = max(alignedShiftMS);
    alignedMeanShiftMS = mean(alignedShiftMS);
    alignedStdShiftMS = std(alignedShiftMS);
    nAlignedClipped = sum(erpAlignedOut.wasClipped);

    fprintf('\nERPAligned shift summary for %s:\n',electrodeName);
    fprintf('Min shift  = %d samples, %.4f ms\n', ...
        min(erpAlignedOut.shiftSamples),alignedMinShiftMS);
    fprintf('Max shift  = %d samples, %.4f ms\n', ...
        max(erpAlignedOut.shiftSamples),alignedMaxShiftMS);
    fprintf('Mean shift = %.2f samples, %.4f ms\n', ...
        mean(erpAlignedOut.shiftSamples),alignedMeanShiftMS);
    fprintf('Std shift  = %.2f samples, %.4f ms\n', ...
        std(erpAlignedOut.shiftSamples),alignedStdShiftMS);
    fprintf('Clipped trials = %d out of %d\n', ...
        nAlignedClipped,length(erpAlignedOut.wasClipped));

    summaryRows(end+1,:) = {electrodeName, ...
        rawError,erpError,shiftError,scaleError,shiftScaleError,alignedError, ...
        erpImprovement,shiftImprovement,scaleImprovement,shiftScaleImprovement,alignedImprovement, ...
        alignedMinShiftMS,alignedMaxShiftMS,alignedMeanShiftMS,alignedStdShiftMS,nAlignedClipped, ...
        erpOutsideChange,shiftedOutsideChange,scaledOutsideChange,shiftScaleOutsideChange,alignedOutsideChange};

    %% Figure 1: ERPAligned time-domain diagnostic
    numTrialsToShow = 30;
    trialSubset = 1:min(numTrialsToShow,size(dataHighStim,1));

    figure('Color','w','Position',[100 100 1100 950]);

    subplot(5,1,1);
    plot(timeVals,dataHighStim(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(dataHighStim,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['Uncleaned high-stimulation trials: ' electrodeName]);
    box off;

    subplot(5,1,2);
    plot(timeVals,erpAlignedOut.alignedTrials(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(erpAlignedOut.alignedTrials,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title('Trials after alignment');
    box off;

    subplot(5,1,3);
    plot(timeVals,erpAlignedOut.template,'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('ERP template computed after trial alignment');
    box off;

    subplot(5,1,4);
    plot(timeVals,erpAlignedOut.trialTemplates(trialSubset,:)','Color',[0.70 0.70 0.70]);
    hold on;
    plot(timeVals,mean(erpAlignedOut.trialTemplates,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('Trial-specific templates shifted back to original timing');
    box off;

    subplot(5,1,5);
    plot(timeVals,erpAlignedOut.cleanedData(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(erpAlignedOut.cleanedData,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['ERPAligned residual trials: ' electrodeName]);
    box off;

    sgtitle(['ERPAligned template subtraction: ' electrodeName]);

    outFileTime = fullfile(figFolder,['demo05_aligned_time_' electrodeName '.png']);
    exportgraphics(gcf,outFileTime,'Resolution',300);

    %% Figure 2: FFT comparison across ERP methods
    figure('Color','w','Position',[100 100 1000 600]);

    plot(freqVals,fftNoStim,'k','LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'r','LineWidth',1.7);
    plot(freqVals,fftERP,'b','LineWidth',1.5);
    plot(freqVals,fftShift,'m','LineWidth',1.5);
    plot(freqVals,fftScale,'c','LineWidth',1.5);
    plot(freqVals,fftShiftScale,'g','LineWidth',1.5);
    plot(freqVals,fftAligned,'LineWidth',1.8);

    xlim([0 200]);

    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');
    title(['FFT comparison across ERP-based methods: ' electrodeName]);

    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERPSubtraction', ...
            'ERPShifted', ...
            'ERPScaled', ...
            'ERPShiftedScaled', ...
            'ERPAligned'}, ...
            'Location','best');

    box off;

    outFileFFT = fullfile(figFolder,['demo05_erp_aligned_fft_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileFFT,'ContentType','vector');

    %% Figure 3: ERPAligned shift diagnostics
    figure('Color','w','Position',[100 100 850 650]);

    subplot(2,1,1);
    histogram(alignedShiftMS,'BinWidth',0.5);
    hold on;
    xline(0,'k--','LineWidth',1.2);
    xline(alignedMeanShiftMS,'r--','LineWidth',1.2);
    title(['ERPAligned shift distribution: ' electrodeName]);
    xlabel('Estimated alignment shift (ms)');
    ylabel('Number of trials');
    legend({'Shift distribution','Zero shift','Mean shift'},'Location','best');
    box off;

    subplot(2,1,2);
    histogram(erpAlignedOut.peakTimes * 1000,'BinWidth',0.5);
    hold on;
    xline(erpAlignedOut.targetPeakTime * 1000,'r--','LineWidth',1.2);
    title('Detected first-pulse timing before alignment');
    xlabel('Detected first-pulse time (ms)');
    ylabel('Number of trials');
    legend({'Detected peaks','Target alignment time'},'Location','best');
    box off;

    outFileDiag = fullfile(figFolder,['demo05_aligned_diagnostics_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileDiag,'ContentType','vector');

    fprintf('\nSaved:\n');
    fprintf('%s\n',outFileTime);
    fprintf('%s\n',outFileFFT);
    fprintf('%s\n',outFileDiag);

end

%% Save summary table
summaryTable = cell2table(summaryRows, ...
    'VariableNames',{'Electrode', ...
    'RawFFTError','ERPFFTError','ShiftedFFTError','ScaledFFTError','ShiftedScaledFFTError','AlignedFFTError', ...
    'ERPImprovementPercent','ShiftedImprovementPercent','ScaledImprovementPercent','ShiftedScaledImprovementPercent','AlignedImprovementPercent', ...
    'AlignedMinShiftMS','AlignedMaxShiftMS','AlignedMeanShiftMS','AlignedStdShiftMS','AlignedNumClippedTrials', ...
    'ERPOutsideWindowChange','ShiftedOutsideWindowChange','ScaledOutsideWindowChange','ShiftedScaledOutsideWindowChange','AlignedOutsideWindowChange'});

disp(summaryTable);

summaryFile = fullfile(figFolder,'demo05_aligned_summary.csv');
writetable(summaryTable,summaryFile);

fprintf('\nSaved ERPAligned summary:\n%s\n',summaryFile);
fprintf('\ndemo_05_ERPAligned complete.\n');