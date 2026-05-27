%% demo_01_ERPSubtraction
% ERP subtraction baseline
% demo for elec1 and elec7

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
plotWindow = [-0.1 0.55];

pos = find(timeVals > artifactWindow(1) & timeVals < artifactWindow(2));

Fs = 1/median(diff(timeVals));
N = length(pos);
freqVals = (0:N-1) * (Fs/N);

freqRangeForMetric = [0 200];
freqMask = freqVals >= freqRangeForMetric(1) & freqVals <= freqRangeForMetric(2);

fprintf('\nERP subtraction demo\n');
fprintf('Artifact/FFT window: %.2f to %.2f s\n',artifactWindow(1),artifactWindow(2));
fprintf('FFT samples: %d\n',N);
fprintf('Frequency resolution: %.2f Hz\n',Fs/N);

%% Electrodes to run
electrodesToRun = {'elec1','elec7'};

summaryRows = {};

for iElec = 1:length(electrodesToRun)

    electrodeName = electrodesToRun{iElec};

    fprintf('\n============================================================\n');
    fprintf('Running ERP subtraction for %s\n',electrodeName);
    fprintf('============================================================\n');

    %% Load electrode data
    lfpFile = fullfile(lfpFolder,[electrodeName '.mat']);
    D = load(lfpFile);

    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    %% Run ERP subtraction
    params = struct();
    params.subtractWindow = artifactWindow;
    params.doBaselineCorrection = false;

    erpOut = ERPSubtraction(dataHighStim,timeVals,params);

    %% FFT calculation
    fftNoStim = log10(mean(abs(fft(dataNoStim(:,pos)'))'));
    fftHighRaw = log10(mean(abs(fft(dataHighStim(:,pos)'))'));
    fftHighERP = log10(mean(abs(fft(erpOut.cleanedData(:,pos)'))'));

    %% Quantitative FFT error
    rawError = norm(fftHighRaw(freqMask) - fftNoStim(freqMask));
    erpError = norm(fftHighERP(freqMask) - fftNoStim(freqMask));

    improvementPercent = 100 * (rawError - erpError) / rawError;

    fprintf('Raw FFT error vs no-stim      = %.4f\n',rawError);
    fprintf('ERP FFT error vs no-stim      = %.4f\n',erpError);
    fprintf('Improvement                   = %.2f %%\n',improvementPercent);

    summaryRows(end+1,:) = {electrodeName,rawError,erpError,improvementPercent};

    %% Time-domain ERP diagnostic figure
    numTrialsToShow = 30;
    trialSubset = 1:min(numTrialsToShow,size(dataHighStim,1));

    meanHighStim = mean(dataHighStim,1);
    meanCleanERP = mean(erpOut.cleanedData,1);

    figure('Color','w','Position',[100 100 1100 750]);

    subplot(3,1,1);
    plot(timeVals,dataHighStim(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,meanHighStim,'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['Uncleaned high-stimulation trials: ' electrodeName]);
    box off;

    subplot(3,1,2);
    plot(timeVals,erpOut.template,'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('Template amplitude (a.u.)');
    title('ERP/template component subtracted');
    box off;

    subplot(3,1,3);
    plot(timeVals,erpOut.cleanedData(trialSubset,:)','Color',[0.75 0.75 0.75]);
    hold on;
    plot(timeVals,meanCleanERP,'k','LineWidth',1.8);
    xlim(plotWindow);
    xlabel('Time (s)');
    ylabel('LFP amplitude (a.u.)');
    title(['ERP-subtracted high-stimulation trials: ' electrodeName]);
    box off;

    sgtitle(['ERP-based template subtraction: ' electrodeName]);

    outFileTime = fullfile(figFolder,['demo01_erp_time_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileTime,'ContentType','vector');

    %% FFT comparison figure
    figure('Color','w','Position',[100 100 900 550]);

    plot(freqVals,fftNoStim,'LineWidth',1.7);
    hold on;
    plot(freqVals,fftHighRaw,'LineWidth',1.7);
    plot(freqVals,fftHighERP,'LineWidth',1.7);

    xlim([0 200]);

    xlabel('Frequency (Hz)');
    ylabel('log_{10} mean FFT magnitude (a.u.)');
    title(['FFT comparison after ERP subtraction: ' electrodeName]);

    legend({'No-stimulation reference', ...
            'Uncleaned 64 \muA', ...
            'ERP-subtracted 64 \muA'}, ...
            'Location','best');

    box off;

    outFileFFT = fullfile(figFolder,['demo01_erp_fft_' electrodeName '.pdf']);
    exportgraphics(gcf,outFileFFT,'ContentType','vector');

    fprintf('Saved:\n');
    fprintf('%s\n',outFileTime);
    fprintf('%s\n',outFileFFT);

end

%% Save summary table
summaryTable = cell2table(summaryRows, ...
    'VariableNames',{'Electrode','RawFFTError','ERPFFTError','ImprovementPercent'});

disp(summaryTable);

summaryFile = fullfile(figFolder,'demo01_erp_summary.csv');
writetable(summaryTable,summaryFile);

fprintf('\nSaved ERP summary:\n%s\n',summaryFile);