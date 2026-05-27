%% demo_04_ERPShiftedScaled
% Compare ERP-based variants:
%   1. No-stimulation reference
%   2. Uncleaned high-stimulation signal
%   3. ERPSubtraction
%   4. ERPShifted
%   5. ERPScaled
%   6. ERPShiftedScaled
%
% Runs elec1 and elec7.
% Uses final window convention:
%   [0 0.4] means timeVals > 0 and timeVals <= 0.4
%   Actual FFT window: 0.0005 to 0.4000 s
%   FFT samples: 800

clear; close all; clc;

clear ERPSubtraction ERPShifted ERPScaled ERPShiftedScaled

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
fprintf('demo_04_ERPShiftedScaled\n');
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
    fprintf('Running ERPShiftedScaled comparison for %s\n',electrodeName);
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

    fprintf('\nSanity checks:\n');
    fprintf('ERPSubtraction max outside-window change    = %.12f\n',erpOutsideChange);
    fprintf('ERPShifted max outside-window change        = %.12f\n',shiftedOutsideChange);
    fprintf('ERPScaled max outside-window change         = %.12f\n',scaledOutsideChange);
    fprintf('ERPShiftedScaled max outside-window change  = %.12f\n',shiftScaleOutsideChange);

    %% FFT calculation
    fftNoStim = log10(mean(abs(fft(dataNoStim(:,pos)'))'));
    fftHighRaw = log10(mean(abs(fft(dataHighStim(:,pos)'))'));
    fftERP = log10(mean(abs(fft(erpOut.cleanedData(:,pos)'))'));
    fftShift = log10(mean(abs(fft(erpShiftOut.cleanedData(:,pos)'))'));
    fftScale = log10(mean(abs(fft(erpScaleOut.cleanedData(:,pos)'))'));
    fftShiftScale = log10(mean(abs(fft(erpShiftScaleOut.cleanedData(:,pos)'))'));

    %% Quantitative FFT errors
    rawError = norm(fftHighRaw(freqMask) - fftNoStim(freqMask));
    erpError = norm(fftERP(freqMask) - fftNoStim(freqMask));
    shiftError = norm(fftShift(freqMask) - fftNoStim(freqMask));
    scaleError = norm(fftScale(freqMask) - fftNoStim(freqMask));
    shiftScaleError = norm(fftShiftScale(freqMask) - fftNoStim(freqMask));

    erpImprovement = 100 * (rawError - erpError) / rawError;
    shiftImprovement = 100 * (rawError - shiftError) / rawError;
    scaleImprovement = 100 * (rawError - scaleError) / rawError;
    shiftScaleImprovement = 100 * (rawError - shiftScaleError) / rawError;

    fprintf('\nFFT comparison results for %s:\n',electrodeName);
    fprintf('Raw FFT error vs no-stim         = %.4f\n',rawError);
    fprintf('ERPSubtraction error             = %.4f\n',erpError);
    fprintf('ERPShifted error                 = %.4f\n',shiftError);
    fprintf('ERPScaled error                  = %.4f\n',scaleError);
    fprintf('ERPShiftedScaled error           = %.4f\n',shiftScaleError);

    fprintf('\nImprovement over uncleaned high-stim:\n');
    fprintf('ERPSubtraction improvement       = %.2f %%\n',erpImprovement);
    fprintf('ERPShifted improvement           = %.2f %%\n',shiftImprovement);
    fprintf('ERPScaled improvement            = %.2f %%\n',scaleImprovement);
    fprintf('ERPShiftedScaled improvement     = %.2f %%\n',shiftScaleImprovement);

    %% Shift and scale summaries
    shiftMS = erpShiftScaleOut.shiftTimes * 1000;
    scaleFactors = erpShiftScaleOut.scaleFactors;

    maxShiftSamples = round((shiftScaleParams.maxShiftMS/1000) / median(diff(timeVals)));

    minShiftMS = min(shiftMS);
    maxShiftMS = max(shiftMS);
    meanShiftMS = mean(shiftMS);
    stdShiftMS = std(shiftMS);

    minScale = min(scaleFactors);
    maxScale = max(scaleFactors);
    meanScale = mean(scaleFactors);
    stdScale = std(scaleFactors);

    nShiftLower = sum(erpShiftScaleOut.shiftSamples == -maxShiftSamples);
    nShiftUpper = sum(erpShiftScaleOut.shiftSamples == maxShiftSamples);

    nScaleLower = sum(scaleFactors == shiftScaleParams.scaleBounds(1));
    nScaleUpper = sum(scaleFactors == shiftScaleParams.scaleBounds(2));

    fprintf('\nShift summary for ERPShiftedScaled, %s:\n',electrodeName);
    fprintf('Min shift  = %d samples, %.4f ms\n', ...
        min(erpShiftScaleOut.shiftSamples),minShiftMS);
    fprintf('Max shift  = %d samples, %.4f ms\n', ...
        max(erpShiftScaleOut.shiftSamples),maxShiftMS);
    fprintf('Mean shift = %.2f samples, %.4f ms\n', ...
        mean(erpShiftScaleOut.shiftSamples),meanShiftMS);
    fprintf('Std shift  = %.2f samples, %.4f ms\n', ...
        std(erpShiftScaleOut.shiftSamples),stdShiftMS);
    fprintf('Number of shifts at lower bound = %d\n',nShiftLower);
    fprintf('Number of shifts at upper bound = %d\n',nShiftUpper);

    fprintf('\nScale summary for ERPShiftedScaled, %s:\n',electrodeName);
    fprintf('Min scale  = %.4f\n',minScale);
    fprintf('Max scale  = %.4f\n',maxScale);
    fprintf('Mean scale = %.4f\n',meanScale);
    fprintf('Std scale  = %.4f\n',stdScale);
    fprintf('Number of scale factors at lower bound %.2f = %d\n', ...
        shiftScaleParams.scaleBounds(1),nScaleLower);
    fprintf('Number of scale factors at upper bound %.2f = %d\n', ...
        shiftScaleParams.scaleBounds(2),nScaleUpper);

    summaryRows(end+1,:) = {electrodeName, ...
        rawError,erpError,shiftError,scaleError,shiftScaleError, ...
        erpImprovement,shiftImprovement,scaleImprovement,shiftScaleImprovement, ...
        minShiftMS,maxShiftMS,meanShiftMS,stdShiftMS,nShiftLower,nShiftUpper, ...
        minScale,maxScale,meanScale,stdScale,nScaleLower,nScaleUpper, ...
        erpOutsideChange,shiftedOutsideChange,scaledOutsideChange,shiftScaleOutsideChange};

    %% Figure 1: Time-domain diagnostic for shifted + scaled
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
    plot(timeVals,erpShiftScaleOut.template,'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('ERP template');
    box off;

    subplot(4,1,3);
    plot(timeVals,erpShiftScaleOut.shiftedScaledTemplates(trialSubset,:)','Color',[0.70 0.70 0.70]);
    hold on;
    plot(timeVals,mean(erpShiftScaleOut.shiftedScaledTemplates,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('Trial-wise shifted + scaled ERP templates');
    box off;

    subplot(4,1,4);
    plot(timeVals,erpShiftScaleOut.cleanedData(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(erpShiftScaleOut.cleanedData,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['ERPShiftedScaled residual trials: ' electrodeName]);
    box off;

    sgtitle(['ERPShiftedScaled template subtraction: ' electrodeName]);

    outFileTime = fullfile(figFolder,['demo04_shifted_scaled_time_' electrodeName '.png']);
    exportgraphics(gcf,outFileTime,'Resolution',300);

    %% Figure 2: FFT comparison across ERP methods
    figure('Color','w','Position',[100 100 1000 600]);

    plot(freqVals,fftNoStim,'k','LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'r','LineWidth',1.7);
    plot(freqVals,fftERP,'b','LineWidth',1.5);
    plot(freqVals,fftShift,'m','LineWidth',1.5);
    plot(freqVals,fftScale,'c','LineWidth',1.5);
    plot(freqVals,fftShiftScale,'g','LineWidth',1.7);

    xlim([0 200]);

    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');
    title(['FFT comparison across ERP-based methods: ' electrodeName]);

    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERPSubtraction', ...
            'ERPShifted', ...
            'ERPScaled', ...
            'ERPShiftedScaled'}, ...
            'Location','best');

    box off;

    outFileFFT = fullfile(figFolder,['demo04_erp_variants_fft_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileFFT,'ContentType','vector');

    %% Figure 3: Shift and scale diagnostics
    figure('Color','w','Position',[100 100 900 750]);

    subplot(3,1,1);
    histogram(shiftMS,'BinWidth',0.5);
    hold on;
    xline(0,'k--','LineWidth',1.2);
    xline(meanShiftMS,'r--','LineWidth',1.2);
    title(['ERPShiftedScaled shift distribution: ' electrodeName]);
    xlabel('Shift (ms)');
    ylabel('Number of trials');
    legend({'Shifts','Zero','Mean'},'Location','best');
    box off;

    subplot(3,1,2);
    histogram(scaleFactors,'BinWidth',0.025);
    hold on;
    xline(1,'k--','LineWidth',1.2);
    xline(meanScale,'r--','LineWidth',1.2);
    title('ERPShiftedScaled scale-factor distribution');
    xlabel('Scale factor');
    ylabel('Number of trials');
    legend({'Scales','Scale = 1','Mean'},'Location','best');
    box off;

    subplot(3,1,3);
    scatter(shiftMS,scaleFactors,25,'filled');
    hold on;
    xline(0,'k--','LineWidth',1.2);
    yline(1,'k--','LineWidth',1.2);
    title('Shift vs scale factor');
    xlabel('Shift (ms)');
    ylabel('Scale factor');
    box off;

    outFileDiag = fullfile(figFolder,['demo04_shifted_scaled_diagnostics_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileDiag,'ContentType','vector');

    fprintf('\nSaved:\n');
    fprintf('%s\n',outFileTime);
    fprintf('%s\n',outFileFFT);
    fprintf('%s\n',outFileDiag);

end

%% Save summary table
summaryTable = cell2table(summaryRows, ...
    'VariableNames',{'Electrode', ...
    'RawFFTError','ERPFFTError','ShiftedFFTError','ScaledFFTError','ShiftedScaledFFTError', ...
    'ERPImprovementPercent','ShiftedImprovementPercent','ScaledImprovementPercent','ShiftedScaledImprovementPercent', ...
    'MinShiftMS','MaxShiftMS','MeanShiftMS','StdShiftMS','NumShiftLowerBound','NumShiftUpperBound', ...
    'MinScale','MaxScale','MeanScale','StdScale','NumScaleLowerBound','NumScaleUpperBound', ...
    'ERPOutsideWindowChange','ShiftedOutsideWindowChange','ScaledOutsideWindowChange','ShiftedScaledOutsideWindowChange'});

disp(summaryTable);

summaryFile = fullfile(figFolder,'demo04_shifted_scaled_summary.csv');
writetable(summaryTable,summaryFile);

fprintf('\nSaved ERPShiftedScaled summary:\n%s\n',summaryFile);
fprintf('\ndemo_04_ERPShiftedScaled complete.\n');