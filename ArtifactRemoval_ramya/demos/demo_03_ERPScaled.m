%% demo_03_ERPScaled
% ERPScaled comparison with ERPSubtraction baseline
% demo for elec1 and elec7
%
% Methods compared:
%   1. No-stimulation reference
%   2. Uncleaned high-stimulation signal
%   3. ERPSubtraction
%   4. ERPScaled

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
fprintf('demo_03_ERPScaled\n');
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
    fprintf('Running ERPScaled comparison for %s\n',electrodeName);
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

    %% Method 2: ERPScaled
    scaleParams = struct();
    scaleParams.subtractWindow = artifactWindow;

    % Scale factor is estimated inside the same artifact/analysis window.
    scaleParams.scaleWindow = artifactWindow;

    % Keep scale factors in a reasonable range.
    scaleParams.scaleBounds = [0 2];

    scaleParams.doBaselineCorrection = false;

    erpScaleOut = ERPScaled(dataHighStim,timeVals,scaleParams);

    %% Sanity checks
    subtractIdx = timeVals > artifactWindow(1) & timeVals <= artifactWindow(2);
    outsideIdx = ~subtractIdx;

    erpOutsideChange = max(abs(erpOut.cleanedData(:,outsideIdx) - ...
                               dataHighStim(:,outsideIdx)),[],'all');

    scaledOutsideChange = max(abs(erpScaleOut.cleanedData(:,outsideIdx) - ...
                                  dataHighStim(:,outsideIdx)),[],'all');

    fprintf('\nSanity checks:\n');
    fprintf('ERPSubtraction max outside-window change = %.12f\n',erpOutsideChange);
    fprintf('ERPScaled max outside-window change      = %.12f\n',scaledOutsideChange);

    %% FFT calculation
    fftNoStim = log10(mean(abs(fft(dataNoStim(:,pos)'))'));
    fftHighRaw = log10(mean(abs(fft(dataHighStim(:,pos)'))'));
    fftERPClean = log10(mean(abs(fft(erpOut.cleanedData(:,pos)'))'));
    fftScaleClean = log10(mean(abs(fft(erpScaleOut.cleanedData(:,pos)'))'));

    %% Quantitative FFT error
    rawError = norm(fftHighRaw(freqMask) - fftNoStim(freqMask));
    erpError = norm(fftERPClean(freqMask) - fftNoStim(freqMask));
    scaleError = norm(fftScaleClean(freqMask) - fftNoStim(freqMask));

    erpImprovement = 100 * (rawError - erpError) / rawError;
    scaleImprovement = 100 * (rawError - scaleError) / rawError;

    fprintf('\nFFT comparison results for %s:\n',electrodeName);
    fprintf('Raw FFT error vs no-stim         = %.4f\n',rawError);
    fprintf('ERPSubtraction error             = %.4f\n',erpError);
    fprintf('ERPScaled error                  = %.4f\n',scaleError);
    fprintf('ERPSubtraction improvement       = %.2f %%\n',erpImprovement);
    fprintf('ERPScaled improvement            = %.2f %%\n',scaleImprovement);

    %% Scale summary
    scaleFactors = erpScaleOut.scaleFactors;

    minScale = min(scaleFactors);
    maxScale = max(scaleFactors);
    meanScale = mean(scaleFactors);
    stdScale = std(scaleFactors);

    nLowerBound = sum(scaleFactors == scaleParams.scaleBounds(1));
    nUpperBound = sum(scaleFactors == scaleParams.scaleBounds(2));

    fprintf('\nScale factor summary for %s:\n',electrodeName);
    fprintf('Min scale  = %.4f\n',minScale);
    fprintf('Max scale  = %.4f\n',maxScale);
    fprintf('Mean scale = %.4f\n',meanScale);
    fprintf('Std scale  = %.4f\n',stdScale);
    fprintf('Number at lower bound %.2f = %d\n',scaleParams.scaleBounds(1),nLowerBound);
    fprintf('Number at upper bound %.2f = %d\n',scaleParams.scaleBounds(2),nUpperBound);

    summaryRows(end+1,:) = {electrodeName, ...
        rawError,erpError,scaleError, ...
        erpImprovement,scaleImprovement, ...
        minScale,maxScale,meanScale,stdScale, ...
        nLowerBound,nUpperBound, ...
        erpOutsideChange,scaledOutsideChange};

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
    plot(timeVals,erpScaleOut.template,'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('ERP template');
    box off;

    subplot(4,1,3);
    plot(timeVals,erpScaleOut.scaledTemplates(trialSubset,:)','Color',[0.70 0.70 0.70]);
    hold on;
    plot(timeVals,mean(erpScaleOut.scaledTemplates,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('Trial-wise scaled ERP templates');
    box off;

    subplot(4,1,4);
    plot(timeVals,erpScaleOut.cleanedData(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(erpScaleOut.cleanedData,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['ERPScaled residual trials: ' electrodeName]);
    box off;

    sgtitle(['ERPScaled template subtraction: ' electrodeName]);

    outFileTime = fullfile(figFolder,['demo03_scaled_time_' electrodeName '.png']);
    exportgraphics(gcf,outFileTime,'Resolution',300);

    %% Figure 2: FFT comparison
    figure('Color','w','Position',[100 100 900 550]);

    plot(freqVals,fftNoStim,'k','LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'r','LineWidth',1.7);
    plot(freqVals,fftERPClean,'b','LineWidth',1.7);
    plot(freqVals,fftScaleClean,'m','LineWidth',1.7);

    xlim([0 200]);

    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');
    title(['FFT comparison: ERPSubtraction vs ERPScaled, ' electrodeName]);

    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERPSubtraction', ...
            'ERPScaled'}, ...
            'Location','best');

    box off;

    outFileFFT = fullfile(figFolder,['demo03_scaled_fft_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileFFT,'ContentType','vector');

    %% Figure 3: Scale distribution
    figure('Color','w','Position',[100 100 850 600]);

    subplot(2,1,1);
    histogram(scaleFactors,'BinWidth',0.025);
    hold on;
    xline(1,'k--','LineWidth',1.2);
    xline(meanScale,'r--','LineWidth',1.2);
    title(['Distribution of estimated ERP scale factors: ' electrodeName]);
    xlabel('Scale factor');
    ylabel('Number of trials');
    legend({'Scale distribution','Scale = 1','Mean scale'},'Location','best');
    box off;

    subplot(2,1,2);
    plot(scaleFactors,'o-','LineWidth',1.0);
    hold on;
    yline(1,'k--','LineWidth',1.2);
    title('Estimated scale factor for each high-stimulation trial');
    xlabel('Trial number');
    ylabel('Scale factor');
    box off;

    outFileScale = fullfile(figFolder,['demo03_scaled_scale_hist_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileScale,'ContentType','vector');

    fprintf('\nSaved:\n');
    fprintf('%s\n',outFileTime);
    fprintf('%s\n',outFileFFT);
    fprintf('%s\n',outFileScale);

end

%% Save summary table
summaryTable = cell2table(summaryRows, ...
    'VariableNames',{'Electrode', ...
    'RawFFTError','ERPFFTError','ScaledFFTError', ...
    'ERPImprovementPercent','ScaledImprovementPercent', ...
    'MinScale','MaxScale','MeanScale','StdScale', ...
    'NumAtLowerBound','NumAtUpperBound', ...
    'ERPOutsideWindowChange','ScaledOutsideWindowChange'});

disp(summaryTable);

summaryFile = fullfile(figFolder,'demo03_scaled_summary.csv');
writetable(summaryTable,summaryFile);

fprintf('\nSaved ERPScaled summary:\n%s\n',summaryFile);
fprintf('\ndemo_03_ERPScaled complete.\n');