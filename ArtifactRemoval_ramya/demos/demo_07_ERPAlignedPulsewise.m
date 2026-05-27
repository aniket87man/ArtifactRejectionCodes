%% demo_07_ERPAlignedPulsewise
%
% This compares:
%   1. No-stimulation reference trials
%   2. Uncleaned high-stimulation trials
%   3. ERP template subtraction
%   4. Aligned ERP template subtraction
%   5. Pulse-wise template subtraction
%   6. ERPAlignedPulsewise hybrid subtraction
%
% The analysis is performed for selected electrodes using the 0--0.4 s
% artifact-analysis window. FFT summaries are computed from trial-wise FFT
% magnitudes and averaged across trials.

clear; close all; clc;

clear ERPSubtraction ERPAligned PulsewiseTemplate ERPAlignedPulsewise

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

% FFT window: 0.0005 to 0.4000 s, 800 samples at 2000 Hz.
pos = find(timeVals > artifactWindow(1) & timeVals <= artifactWindow(2));

Fs = 1/median(diff(timeVals));
N = length(pos);
freqVals = (0:N-1) * (Fs/N);

freqRangeForMetric = [0 200];
freqMask = freqVals >= freqRangeForMetric(1) & freqVals <= freqRangeForMetric(2);

fprintf('\n============================================================\n');
fprintf('demo_07_ERPAlignedPulsewise\n');
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
    fprintf('Running ERPAlignedPulsewise comparison for %s\n',electrodeName);
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

    %% Method 2: ERPAligned
    alignedParams = struct();
    alignedParams.subtractWindow = artifactWindow;
    alignedParams.alignWindow = [-0.01 0.03];
    alignedParams.maxShiftMS = 10;
    alignedParams.doBaselineCorrection = false;

    erpAlignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);

    %% Method 3: PulsewiseTemplate
    pulseParams = struct();

    pulseParams.subtractWindow = artifactWindow;
    pulseParams.pulseSearchWindow = [-0.02 0.32];
    pulseParams.pulseWindowMS = [-5 35];
    pulseParams.minPulseDistanceMS = 25;
    pulseParams.expectedNumPulses = 7;
    pulseParams.maxNumPulses = 20;
    pulseParams.thresholdMAD = 6;
    pulseParams.localBaselineMS = 20;
    pulseParams.templateStatistic = 'median';
    pulseParams.taperEdgeMS = 2;
    pulseParams.doBaselineCorrection = false;

    pulseOut = PulsewiseTemplate(dataHighStim,timeVals,pulseParams);

    %% Method 4: ERPAlignedPulsewise hybrid
    hybridPulseParams = pulseParams;

    % A lower threshold is used because the pulse-wise stage operates on
    % the ERPAligned residual, where pulse residuals can be smaller.
    hybridPulseParams.thresholdMAD = 3;

    hybridParams = struct();
    hybridParams.subtractWindow = artifactWindow;
    hybridParams.alignedParams = alignedParams;
    hybridParams.pulseParams = hybridPulseParams;

    hybridOut = ERPAlignedPulsewise(dataHighStim,timeVals,hybridParams);

    fprintf('\nPulsewiseTemplate detected %d pulses.\n',length(pulseOut.pulseTimes));
    fprintf('PulsewiseTemplate pulse times, ms:\n');
    disp(pulseOut.pulseTimes * 1000);

    fprintf('\nERPAlignedPulsewise residual stage detected %d pulses.\n', ...
        length(hybridOut.pulseOut.pulseTimes));
    fprintf('ERPAlignedPulsewise residual pulse times, ms:\n');
    disp(hybridOut.pulseOut.pulseTimes * 1000);

    %% Sanity checks
    subtractIdx = timeVals > artifactWindow(1) & timeVals <= artifactWindow(2);
    outsideIdx = ~subtractIdx;

    erpOutsideChange = max(abs(erpOut.cleanedData(:,outsideIdx) - ...
                               dataHighStim(:,outsideIdx)),[],'all');

    alignedOutsideChange = max(abs(erpAlignedOut.cleanedData(:,outsideIdx) - ...
                                   dataHighStim(:,outsideIdx)),[],'all');

    pulseOutsideChange = max(abs(pulseOut.cleanedData(:,outsideIdx) - ...
                                 dataHighStim(:,outsideIdx)),[],'all');

    hybridOutsideChange = max(abs(hybridOut.cleanedData(:,outsideIdx) - ...
                                  dataHighStim(:,outsideIdx)),[],'all');

    fprintf('\nSanity checks:\n');
    fprintf('ERPSubtraction max outside-window change       = %.12f\n',erpOutsideChange);
    fprintf('ERPAligned max outside-window change           = %.12f\n',alignedOutsideChange);
    fprintf('PulsewiseTemplate max outside-window change    = %.12f\n',pulseOutsideChange);
    fprintf('ERPAlignedPulsewise max outside-window change  = %.12f\n',hybridOutsideChange);

    %% FFT calculation
    % FFT magnitudes are computed for each trial and then averaged.
    fftNoStim = log10(mean(abs(fft(dataNoStim(:,pos)'))'));
    fftHighRaw = log10(mean(abs(fft(dataHighStim(:,pos)'))'));
    fftERP = log10(mean(abs(fft(erpOut.cleanedData(:,pos)'))'));
    fftAligned = log10(mean(abs(fft(erpAlignedOut.cleanedData(:,pos)'))'));
    fftPulse = log10(mean(abs(fft(pulseOut.cleanedData(:,pos)'))'));
    fftHybrid = log10(mean(abs(fft(hybridOut.cleanedData(:,pos)'))'));

    %% Quantitative FFT errors
    rawError = norm(fftHighRaw(freqMask) - fftNoStim(freqMask));
    erpError = norm(fftERP(freqMask) - fftNoStim(freqMask));
    alignedError = norm(fftAligned(freqMask) - fftNoStim(freqMask));
    pulseError = norm(fftPulse(freqMask) - fftNoStim(freqMask));
    hybridError = norm(fftHybrid(freqMask) - fftNoStim(freqMask));

    erpImprovement = 100 * (rawError - erpError) / rawError;
    alignedImprovement = 100 * (rawError - alignedError) / rawError;
    pulseImprovement = 100 * (rawError - pulseError) / rawError;
    hybridImprovement = 100 * (rawError - hybridError) / rawError;

    fprintf('\nFFT comparison results for %s:\n',electrodeName);
    fprintf('Raw FFT error vs no-stim             = %.4f\n',rawError);
    fprintf('ERPSubtraction error                 = %.4f\n',erpError);
    fprintf('ERPAligned error                     = %.4f\n',alignedError);
    fprintf('PulsewiseTemplate error              = %.4f\n',pulseError);
    fprintf('ERPAlignedPulsewise error            = %.4f\n',hybridError);

    fprintf('\nImprovement over uncleaned high-stim:\n');
    fprintf('ERPSubtraction improvement           = %.2f %%\n',erpImprovement);
    fprintf('ERPAligned improvement               = %.2f %%\n',alignedImprovement);
    fprintf('PulsewiseTemplate improvement        = %.2f %%\n',pulseImprovement);
    fprintf('ERPAlignedPulsewise improvement      = %.2f %%\n',hybridImprovement);

    summaryRows(end+1,:) = {electrodeName, ...
        rawError,erpError,alignedError,pulseError,hybridError, ...
        erpImprovement,alignedImprovement,pulseImprovement,hybridImprovement, ...
        length(pulseOut.pulseTimes),length(hybridOut.pulseOut.pulseTimes), ...
        erpOutsideChange,alignedOutsideChange,pulseOutsideChange,hybridOutsideChange};

    %% Figure 1: Stage-wise hybrid diagnostic
    numTrialsToShow = 30;
    trialSubset = 1:min(numTrialsToShow,size(dataHighStim,1));

    figure('Color','w','Position',[100 100 1100 950]);

    subplot(5,1,1);
    plot(timeVals,dataHighStim(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(dataHighStim,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    title(['Uncleaned high-stimulation trials: ' electrodeName]);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    box off;

    subplot(5,1,2);
    plot(timeVals,hybridOut.erpAlignedOut.cleanedData(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(hybridOut.erpAlignedOut.cleanedData,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    title('After ERPAligned stage');
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    box off;

    subplot(5,1,3);
    plot(timeVals,hybridOut.pulseOut.fullTemplate,'k','LineWidth',1.8);
    xlim(plotWindow);
    title('Pulse-wise template estimated from ERPAligned residual');
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    box off;

    subplot(5,1,4);
    plot(timeVals,hybridOut.cleanedData(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,mean(hybridOut.cleanedData,1),'k','LineWidth',1.8);
    xlim(plotWindow);
    title(['ERPAlignedPulsewise residual trials: ' electrodeName]);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    box off;

    subplot(5,1,5);
    plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.3);
    hold on;
    plot(timeVals,mean(hybridOut.erpAlignedOut.cleanedData,1),'g','LineWidth',1.3);
    plot(timeVals,mean(hybridOut.cleanedData,1),'m','LineWidth',1.3);
    xlim(plotWindow);
    title('Mean signal across processing stages');
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    legend({'Uncleaned 64 \muA','After ERPAligned','After ERPAlignedPulsewise'}, ...
        'Location','best');
    box off;

    sgtitle(['ERPAlignedPulsewise processing stages: ' electrodeName]);

    outFileStages = fullfile(figFolder,['demo07_hybrid_stages_' electrodeName '.png']);
    exportgraphics(gcf,outFileStages,'Resolution',300);

    %% Figure 2: Residual pulse detection
    figure('Color','w','Position',[100 100 1100 750]);

    subplot(3,1,1);
    plot(timeVals,mean(hybridOut.erpAlignedOut.cleanedData,1),'g','LineWidth',1.3);
    hold on;
    plot(hybridOut.pulseOut.pulseTimes, ...
        mean(hybridOut.erpAlignedOut.cleanedData(:,hybridOut.pulseOut.pulseSamples),1), ...
        'ko','MarkerFaceColor','k');
    xlim(plotWindow);
    title(['Residual pulse detection after ERPAligned: ' electrodeName]);
    xlabel('Time (s)');
    ylabel('Residual LFP amplitude (a.u.)');
    box off;

    subplot(3,1,2);
    plot(timeVals,hybridOut.pulseOut.detectionScore,'LineWidth',1.3);
    hold on;
    yline(hybridOut.pulseOut.detectionThreshold,'r--','LineWidth',1.2);
    plot(hybridOut.pulseOut.pulseTimes, ...
        hybridOut.pulseOut.detectionScore(hybridOut.pulseOut.pulseSamples), ...
        'ko','MarkerFaceColor','k');
    xlim(plotWindow);
    title('Residual pulse detection score');
    xlabel('Time (s)');
    ylabel('Detection score');
    legend({'Detection score','Threshold','Detected pulses'},'Location','best');
    box off;

    subplot(3,1,3);
    plot(timeVals,hybridOut.pulseOut.fullTemplate,'k','LineWidth',1.8);
    xlim(plotWindow);
    title('Residual pulse-wise template');
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    box off;

    outFileDetect = fullfile(figFolder,['demo07_hybrid_residual_detection_' electrodeName '.png']);
    exportgraphics(gcf,outFileDetect,'Resolution',300);

    %% Figure 3: Template comparison
    figure('Color','w','Position',[100 100 1100 850]);

    subplot(4,1,1);
    plot(timeVals,erpOut.template,'b','LineWidth',1.5);
    xlim(plotWindow);
    title('ERPSubtraction template');
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    box off;

    subplot(4,1,2);
    plot(timeVals,erpAlignedOut.template,'g','LineWidth',1.5);
    xlim(plotWindow);
    title('ERPAligned template');
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    box off;

    subplot(4,1,3);
    plot(timeVals,pulseOut.fullTemplate,'m','LineWidth',1.5);
    xlim(plotWindow);
    title('PulsewiseTemplate template');
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    box off;

    subplot(4,1,4);
    plot(timeVals,mean(hybridOut.totalTemplatePerTrial,1),'k','LineWidth',1.5);
    xlim(plotWindow);
    title('ERPAlignedPulsewise average total template');
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    box off;

    sgtitle(['Template comparison: ' electrodeName]);

    outFileTemplates = fullfile(figFolder,['demo07_hybrid_template_comparison_' electrodeName '.png']);
    exportgraphics(gcf,outFileTemplates,'Resolution',300);

    %% Figure 4: FFT comparison
    figure('Color','w','Position',[100 100 1000 600]);

    plot(freqVals,fftNoStim,'k','LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'r','LineWidth',1.7);
    plot(freqVals,fftERP,'b','LineWidth',1.5);
    plot(freqVals,fftAligned,'g','LineWidth',1.5);
    plot(freqVals,fftPulse,'c','LineWidth',1.5);
    plot(freqVals,fftHybrid,'m','LineWidth',1.8);

    xlim([0 200]);

    title(['FFT comparison: ERPAlignedPulsewise, ' electrodeName]);
    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');

    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERPSubtraction', ...
            'ERPAligned', ...
            'PulsewiseTemplate', ...
            'ERPAlignedPulsewise'}, ...
            'Location','best');

    box off;

    outFileFFT = fullfile(figFolder,['demo07_hybrid_fft_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileFFT,'ContentType','vector');

    %% Figure 5: Simple comparison figure
    figure('Color','w','Position',[100 100 1100 750]);

    subplot(2,1,1);
    plot(timeVals,mean(dataNoStim,1),'k','LineWidth',1.7);
    hold on;
    plot(timeVals,mean(dataHighStim,1),'r','LineWidth',1.5);
    plot(timeVals,mean(erpOut.cleanedData,1),'b','LineWidth',1.5);
    plot(timeVals,mean(hybridOut.cleanedData,1),'m','LineWidth',1.7);

    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['Time-domain comparison: ERPAlignedPulsewise, ' electrodeName]);
    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERPSubtraction', ...
            'ERPAlignedPulsewise'}, ...
            'Location','best');
    box off;

    subplot(2,1,2);
    plot(freqVals,fftNoStim,'k','LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'r','LineWidth',1.5);
    plot(freqVals,fftERP,'b','LineWidth',1.5);
    plot(freqVals,fftHybrid,'m','LineWidth',1.8);

    xlim([0 200]);
    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');
    title(['Frequency-domain comparison: ERPAlignedPulsewise, ' electrodeName]);
    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERPSubtraction', ...
            'ERPAlignedPulsewise'}, ...
            'Location','best');
    box off;

    sgtitle(['ERPAlignedPulsewise comparison with ERPSubtraction: ' electrodeName]);

    outFileSimple = fullfile(figFolder,['demo07_simple_erp_vs_hybrid_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileSimple,'ContentType','vector');

    fprintf('\nSaved:\n');
    fprintf('%s\n',outFileStages);
    fprintf('%s\n',outFileDetect);
    fprintf('%s\n',outFileTemplates);
    fprintf('%s\n',outFileFFT);
    fprintf('%s\n',outFileSimple);

end

%% Save summary table
summaryTable = cell2table(summaryRows, ...
    'VariableNames',{'Electrode', ...
    'RawFFTError','ERPFFTError','AlignedFFTError','PulsewiseFFTError','HybridFFTError', ...
    'ERPImprovementPercent','AlignedImprovementPercent','PulsewiseImprovementPercent','HybridImprovementPercent', ...
    'NumPulsewisePulses','NumHybridResidualPulses', ...
    'ERPOutsideWindowChange','AlignedOutsideWindowChange','PulsewiseOutsideWindowChange','HybridOutsideWindowChange'});

disp(summaryTable);

summaryFile = fullfile(figFolder,'demo07_hybrid_summary.csv');
writetable(summaryTable,summaryFile);

fprintf('\nSaved ERPAlignedPulsewise summary:\n%s\n',summaryFile);
fprintf('\ndemo_07_ERPAlignedPulsewise complete.\n');