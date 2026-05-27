clear; close all; clc;

%% demo_25_TimeSeriesBeforeAfterCleaning
%
% Purpose:
% Visualize trial-wise LFP time series before and after artifact removal.
%
% Methods shown:
%   1. ERPSubtraction
%   2. PCATemplate K=10
%   3. SMARTALite Ensemble K3K5
%
% Figures per electrode:
%   1. Trial overlay comparison
%   2. Mean before/after comparison
%
% Both full-window and zoomed-window figures are saved.
%
% Interpretation:
%   ERPSubtraction is the conservative baseline.
%   PCATemplate K=10 was strongest in semi-synthetic known-clean validation.
%   SMARTALite Ensemble K3K5 was competitive in real-data total FFT error.

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName  = 'dona';
gridType     = 'Microelectrode';
expDate      = '290825';
protocolName = 'GRF_001';

% Conditions
noStimCondition   = {1,1,1,5,5,4};
highStimCondition = {7,1,1,5,5,4};

% Representative electrodes
% elec4 and elec8 were already selected earlier.
% elec7 is added because it is a useful representative analyzed electrode.
electrodesToPlot = [4 7 8];

% Windows to display
plotWindows = {[-0.10 0.55], [0.02 0.38]};
plotWindowNames = {'full','zoom'};

% Artifact / cleaning window
artifactWindow = [0 0.4];

% SMARTALite settings
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Number of trials to show in trial-overlay panels
numTrialsToShow = 40;

saveFigures = true;

%% ========================================================================
% Paths and loading
% ========================================================================

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder   = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

I = load(lfpInfoFile);
P = load(paramFile);

timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

noStimTrials = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));

%% ========================================================================
% Check required functions
% ========================================================================

requiredMethods = {'ERPSubtraction','PCATemplate','SMARTALite'};

fprintf('\nChecking required functions:\n');

for i = 1:length(requiredMethods)
    fn = requiredMethods{i};
    fnPath = which(fn);

    if isempty(fnPath)
        warning('%s not found on MATLAB path.',fn);
    else
        fprintf('%s -> %s\n',fn,fnPath);
    end
end

%% ========================================================================
% Method options
% ========================================================================

% ERPSubtraction
erpParams = struct();
erpParams.subtractWindow = artifactWindow;
erpParams.doBaselineCorrection = false;

% PCATemplate K=10
pcaParamsK10 = struct();
pcaParamsK10.artifactWindow = artifactWindow;
pcaParamsK10.numComponents = 10;
pcaParamsK10.removeMeanTemplate = true;
pcaParamsK10.taperEdgeMS = 2;
pcaParamsK10.doBaselineCorrection = false;

% SMARTALite K=3
smartaK3Opts = struct();
smartaK3Opts.pulseTimes = pulseTimes;
smartaK3Opts.stiFreq = smartaStimFreq;
smartaK3Opts.prePulse = smartaPrePulse;
smartaK3Opts.K = 3;
smartaK3Opts.window = artifactWindow;
smartaK3Opts.computeFFT = false;

% SMARTALite K=5
smartaK5Opts = smartaK3Opts;
smartaK5Opts.K = 5;

%% ========================================================================
% Output folder
% ========================================================================

figFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','figures','demo25_TimeSeriesBeforeAfterCleaning');

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Main loop
% ========================================================================

for iElec = 1:length(electrodesToPlot)

    elecNum = electrodesToPlot(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    fprintf('\n============================================================\n');
    fprintf('Processing elec%d\n',elecNum);
    fprintf('============================================================\n');

    if ~exist(elecFile,'file')
        warning('Could not find %s. Skipping.',elecFile);
        continue;
    end

    D = load(elecFile);

    if ~isfield(D,'analogData')
        warning('analogData not found in elec%d. Skipping.',elecNum);
        continue;
    end

    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    %% Apply comparison methods

    fprintf('Running ERPSubtraction...\n');
    erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);
    erpCleaned = getCleanedData(erpOut);

    fprintf('Running PCATemplate K=10...\n');
    pcaOut = PCATemplate(dataHighStim,timeVals,pcaParamsK10);
    pcaCleaned = getCleanedData(pcaOut);

    fprintf('Running SMARTALite K=3...\n');
    smartaK3Out = SMARTALite(dataHighStim,timeVals,smartaK3Opts);
    smartaK3Cleaned = getCleanedData(smartaK3Out);

    fprintf('Running SMARTALite K=5...\n');
    smartaK5Out = SMARTALite(dataHighStim,timeVals,smartaK5Opts);
    smartaK5Cleaned = getCleanedData(smartaK5Out);

    ensembleCleaned = 0.5 * (smartaK3Cleaned + smartaK5Cleaned);

    %% Make figures for full and zoomed windows

    for iWin = 1:length(plotWindows)

        plotWindow = plotWindows{iWin};
        windowName = plotWindowNames{iWin};

        plotIdx = timeVals >= plotWindow(1) & timeVals <= plotWindow(2);

        % Use cleaned/no-stim scale so the cleaned traces remain visible.
        % Raw high-stim may clip, which is acceptable for showing artifact strength.
        scaleData = [ ...
            dataNoStim(:,plotIdx); ...
            erpCleaned(:,plotIdx); ...
            pcaCleaned(:,plotIdx); ...
            ensembleCleaned(:,plotIdx)];

        yLimVal = prctile(abs(scaleData(:)),99.5);
        yLimVal = max(yLimVal,eps);
        yLims = [-yLimVal yLimVal];

        %% Figure 1: trial-wise before/after comparison

        fig1 = figure('Name',['demo25 elec' num2str(elecNum) ' ' windowName ' trial overlay'], ...
            'Color','w','Position',[100 100 1200 900]);

        tiledlayout(5,1,'TileSpacing','compact','Padding','compact');

        nexttile;
        plotTrialOverlay(timeVals(plotIdx),dataNoStim(:,plotIdx),[0.75 0.75 0.75],numTrialsToShow);
        hold on;
        plot(timeVals(plotIdx),mean(dataNoStim(:,plotIdx),1),'k','LineWidth',2);
        addPulseLines(pulseTimes,yLims,plotWindow);
        ylim(yLims);
        title(['elec' num2str(elecNum) ': no-stimulation trials']);
        ylabel('LFP amplitude (a.u.)');
        grid on;

        nexttile;
        plotTrialOverlay(timeVals(plotIdx),dataHighStim(:,plotIdx),[0.80 0.80 0.80],numTrialsToShow);
        hold on;
        plot(timeVals(plotIdx),mean(dataHighStim(:,plotIdx),1),'Color',[0.35 0.35 0.35],'LineWidth',2);
        addPulseLines(pulseTimes,yLims,plotWindow);
        ylim(yLims);
        title('Raw high-stimulation trials before cleaning');
        ylabel('LFP amplitude (a.u.)');
        grid on;

        nexttile;
        plotTrialOverlay(timeVals(plotIdx),pcaCleaned(:,plotIdx),[0.65 0.75 0.95],numTrialsToShow);
        hold on;
        plot(timeVals(plotIdx),mean(pcaCleaned(:,plotIdx),1),'Color',[0 0.20 0.80],'LineWidth',2);
        addPulseLines(pulseTimes,yLims,plotWindow);
        ylim(yLims);
        title('Cleaned high-stimulation trials: PCATemplate K=10');
        ylabel('LFP amplitude (a.u.)');
        grid on;

        nexttile;
        plotTrialOverlay(timeVals(plotIdx),ensembleCleaned(:,plotIdx),[0.65 0.85 0.65],numTrialsToShow);
        hold on;
        plot(timeVals(plotIdx),mean(ensembleCleaned(:,plotIdx),1),'Color',[0 0.45 0],'LineWidth',2);
        addPulseLines(pulseTimes,yLims,plotWindow);
        ylim(yLims);
        title('Cleaned high-stimulation trials: SMARTALite ensemble K3K5');
        ylabel('LFP amplitude (a.u.)');
        grid on;

        nexttile;
        plot(timeVals(plotIdx),mean(dataNoStim(:,plotIdx),1),'k','LineWidth',2);
        hold on;
        plot(timeVals(plotIdx),mean(dataHighStim(:,plotIdx),1),'Color',[0.65 0.65 0.65],'LineWidth',1.4);
        plot(timeVals(plotIdx),mean(erpCleaned(:,plotIdx),1),'LineWidth',1.4);
        plot(timeVals(plotIdx),mean(pcaCleaned(:,plotIdx),1),'Color',[0 0.20 0.80],'LineWidth',1.6);
        plot(timeVals(plotIdx),mean(ensembleCleaned(:,plotIdx),1),'Color',[0 0.45 0],'LineWidth',1.8);
        addPulseLines(pulseTimes,yLims,plotWindow);
        ylim(yLims);
        title('Mean comparison');
        xlabel('Time (s)');
        ylabel('Mean LFP amplitude (a.u.)');
        legend({'No-stim','Raw high-stim','ERPSubtraction','PCATemplate K=10','SMARTALite ensemble K3K5'}, ...
            'Location','best');
        grid on;

        sgtitle(['Time-domain artifact-removal comparison: elec' num2str(elecNum) ...
            ', ' windowName ' window']);

        if saveFigures
            baseName1 = ['demo25_elec' num2str(elecNum) '_' windowName '_trial_overlay_methods'];
            saveFigureAllFormats(fig1,figFolder,baseName1);
        end

        %% Figure 2: zoom/full mean-only comparison

        fig2 = figure('Name',['demo25 elec' num2str(elecNum) ' ' windowName ' mean comparison'], ...
            'Color','w','Position',[100 100 1000 550]);

        plot(timeVals(plotIdx),mean(dataNoStim(:,plotIdx),1),'k','LineWidth',2);
        hold on;
        plot(timeVals(plotIdx),mean(dataHighStim(:,plotIdx),1),'Color',[0.65 0.65 0.65],'LineWidth',1.4);
        plot(timeVals(plotIdx),mean(erpCleaned(:,plotIdx),1),'LineWidth',1.4);
        plot(timeVals(plotIdx),mean(pcaCleaned(:,plotIdx),1),'Color',[0 0.20 0.80],'LineWidth',1.6);
        plot(timeVals(plotIdx),mean(ensembleCleaned(:,plotIdx),1),'Color',[0 0.45 0],'LineWidth',1.8);

        addPulseLines(pulseTimes,yLims,plotWindow);

        xlim(plotWindow);
        ylim(yLims);

        xlabel('Time (s)');
        ylabel('Mean LFP amplitude (a.u.)');
        title(['elec' num2str(elecNum) ': mean signal before and after cleaning']);
        legend({'No-stim','Raw high-stim','ERPSubtraction','PCATemplate K=10','SMARTALite ensemble K3K5'}, ...
            'Location','best');
        grid on;

        if saveFigures
            baseName2 = ['demo25_elec' num2str(elecNum) '_' windowName '_mean_comparison_methods'];
            saveFigureAllFormats(fig2,figFolder,baseName2);
        end
    end
end

fprintf('\nSaved demo_25 figures in:\n%s\n',figFolder);

%% ========================================================================
% Local helper functions
% ========================================================================

function cleanedData = getCleanedData(methodOut)

    if isfield(methodOut,'cleanedData')
        cleanedData = methodOut.cleanedData;
    elseif isfield(methodOut,'cleanedTrials')
        cleanedData = methodOut.cleanedTrials;
    else
        error('Method output has neither cleanedData nor cleanedTrials.');
    end
end

function plotTrialOverlay(t,data,colorVal,numTrialsToShow)

    nTrials = size(data,1);
    nShow = min(numTrialsToShow,nTrials);

    if nShow < nTrials
        trialIdx = round(linspace(1,nTrials,nShow));
    else
        trialIdx = 1:nTrials;
    end

    hold on;

    for i = 1:length(trialIdx)
        iTrial = trialIdx(i);
        plot(t,data(iTrial,:),'Color',colorVal,'LineWidth',0.4);
    end
end

function addPulseLines(pulseTimes,yLims,plotWindow)

    for iPulse = 1:length(pulseTimes)

        thisPulse = pulseTimes(iPulse);

        if thisPulse >= plotWindow(1) && thisPulse <= plotWindow(2)
            xline(thisPulse,'--','Color',[0.8 0 0],'LineWidth',0.8);
        end
    end

    ylim(yLims);
end

function saveFigureAllFormats(figHandle,figFolder,baseName)

    pngFile = fullfile(figFolder,[baseName '.png']);
    figFile = fullfile(figFolder,[baseName '.fig']);
    pdfFile = fullfile(figFolder,[baseName '.pdf']);

    saveas(figHandle,pngFile);
    savefig(figHandle,figFile);
    exportgraphics(figHandle,pdfFile,'ContentType','vector');

    fprintf('Saved figure:\n%s\n%s\n%s\n',pngFile,figFile,pdfFile);
end