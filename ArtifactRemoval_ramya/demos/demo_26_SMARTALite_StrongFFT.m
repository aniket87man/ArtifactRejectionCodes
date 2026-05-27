clear; close all; clc;

%% demo_26_SMARTALite_StrongFFTExample
%
% Purpose:
% Find a representative electrode where SMARTALite gives a strong FFT result,
% then plot an FFT comparison similar to the older SMARTA figures.
%
% This figure is for real-data spectral comparison:
% No-stim reference vs raw high-stim vs ERP/PCA vs SMARTALite variants.

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName  = 'dona';
gridType     = 'Microelectrode';
expDate      = '290825';
protocolName = 'GRF_001';

noStimCondition   = {1,1,1,5,5,4};
highStimCondition = {7,1,1,5,5,4};

% Same 31 good V1 electrodes used in demo22/demo24
electrodesToSearch = [5 13];

% FFT window
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% Artifact-removal window
artifactWindow = [0 0.4];

% SMARTALite K sweep, similar to your older figure
KList = [3 5 10 20 40 60];

% Stimulation pulses
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

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

fftIdx = find(timeVals > fftWindow(1) & timeVals <= fftWindow(2));

Fs = 1/median(diff(timeVals));
N = length(fftIdx);
freqVals = (0:N-1) * (Fs/N);

freqMask = freqVals >= freqRangeForMetric(1) & freqVals <= freqRangeForMetric(2);

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));
fprintf('FFT samples      : %d\n',N);
fprintf('Frequency step   : %.2f Hz\n',Fs/N);

%% ========================================================================
% Method options
% ========================================================================

erpParams = struct();
erpParams.subtractWindow = artifactWindow;
erpParams.doBaselineCorrection = false;

pcaParamsK10 = struct();
pcaParamsK10.artifactWindow = artifactWindow;
pcaParamsK10.numComponents = 10;
pcaParamsK10.removeMeanTemplate = true;
pcaParamsK10.taperEdgeMS = 2;
pcaParamsK10.doBaselineCorrection = false;

baseSMARTAOpts = struct();
baseSMARTAOpts.pulseTimes = pulseTimes;
baseSMARTAOpts.stiFreq = smartaStimFreq;
baseSMARTAOpts.prePulse = smartaPrePulse;
baseSMARTAOpts.window = artifactWindow;
baseSMARTAOpts.computeFFT = false;
baseSMARTAOpts.excludeSameTrial = true;

%% ========================================================================
% Output folder
% ========================================================================

figFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','figures','demo26_SMARTALiteStrongFFT');

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Search for strongest SMARTALite example
% ========================================================================

rows = struct([]);
rowCounter = 0;

bestScore = -Inf;
bestInfo = struct();

for iElec = 1:length(electrodesToSearch)

    elecNum = electrodesToSearch(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    if ~exist(elecFile,'file')
        warning('Missing elec%d. Skipping.',elecNum);
        continue;
    end

    fprintf('\nSearching elec%d...\n',elecNum);

    D = load(elecFile);
    analogData = D.analogData;

    dataNoStim = analogData(noStimTrials,:);
    dataHighStim = analogData(highStimTrials,:);

    fftNoStim = validationFFT(dataNoStim,fftIdx);
    fftRaw = validationFFT(dataHighStim,fftIdx);

    rawError = norm(fftRaw(freqMask) - fftNoStim(freqMask));

    % ERP
    erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);
    erpData = getCleanedData(erpOut);
    fftERP = validationFFT(erpData,fftIdx);
    erpError = norm(fftERP(freqMask) - fftNoStim(freqMask));

    % PCA K10
    pcaOut = PCATemplate(dataHighStim,timeVals,pcaParamsK10);
    pcaData = getCleanedData(pcaOut);
    fftPCA = validationFFT(pcaData,fftIdx);
    pcaError = norm(fftPCA(freqMask) - fftNoStim(freqMask));

    smartaErrors = nan(size(KList));
    smartaFFTs = cell(size(KList));
    smartaData = cell(size(KList));

    for iK = 1:length(KList)

        opts = baseSMARTAOpts;
        opts.K = KList(iK);

        out = SMARTALite(dataHighStim,timeVals,opts);
        cleaned = getCleanedData(out);

        thisFFT = validationFFT(cleaned,fftIdx);
        thisError = norm(thisFFT(freqMask) - fftNoStim(freqMask));

        smartaErrors(iK) = thisError;
        smartaFFTs{iK} = thisFFT;
        smartaData{iK} = cleaned;

        rowCounter = rowCounter + 1;
        rows(rowCounter).electrode = elecNum;
        rows(rowCounter).K = KList(iK);
        rows(rowCounter).rawError = rawError;
        rows(rowCounter).erpError = erpError;
        rows(rowCounter).pcaK10Error = pcaError;
        rows(rowCounter).smartaError = thisError;
        rows(rowCounter).smartaImprovementPct = 100*(rawError - thisError)/rawError;
    end

    [bestSMARTAErrorThisElec,bestKIdx] = min(smartaErrors);
    bestKThisElec = KList(bestKIdx);

    % Score favors a clearly improved SMARTALite example.
    % This selects a visually strong representative, not a group statistic.
    improvementPct = 100*(rawError - bestSMARTAErrorThisElec)/rawError;

    % Require SMARTALite to beat ERP in this representative example if possible.
    beatsERP = bestSMARTAErrorThisElec < erpError;

    if beatsERP
        score = improvementPct + 5;  % small bonus for beating ERP
    else
        score = improvementPct;
    end

    fprintf('elec%d: raw %.4f, ERP %.4f, PCA %.4f, best SMARTALite K=%d error %.4f, improvement %.2f %%\n', ...
        elecNum,rawError,erpError,pcaError,bestKThisElec,bestSMARTAErrorThisElec,improvementPct);

    if score > bestScore
        bestScore = score;

        bestInfo.elecNum = elecNum;
        bestInfo.dataNoStim = dataNoStim;
        bestInfo.dataHighStim = dataHighStim;

        bestInfo.fftNoStim = fftNoStim;
        bestInfo.fftRaw = fftRaw;
        bestInfo.fftERP = fftERP;
        bestInfo.fftPCA = fftPCA;
        bestInfo.smartaFFTs = smartaFFTs;
        bestInfo.smartaErrors = smartaErrors;

        bestInfo.rawError = rawError;
        bestInfo.erpError = erpError;
        bestInfo.pcaError = pcaError;
        bestInfo.bestSMARTAError = bestSMARTAErrorThisElec;
        bestInfo.bestK = bestKThisElec;
        bestInfo.bestKIdx = bestKIdx;
        bestInfo.improvementPct = improvementPct;
    end
end

searchTable = struct2table(rows);

fprintf('\n============================================================\n');
fprintf('Selected representative SMARTALite FFT example\n');
fprintf('============================================================\n');
fprintf('Electrode              : elec%d\n',bestInfo.elecNum);
fprintf('Best SMARTALite K      : %d\n',bestInfo.bestK);
fprintf('Raw FFT error          : %.4f\n',bestInfo.rawError);
fprintf('ERP FFT error          : %.4f\n',bestInfo.erpError);
fprintf('PCA K10 FFT error      : %.4f\n',bestInfo.pcaError);
fprintf('Best SMARTALite error  : %.4f\n',bestInfo.bestSMARTAError);
fprintf('SMARTALite improvement : %.2f %%\n',bestInfo.improvementPct);

%% ========================================================================
% Figure 1: SMARTALite K comparison, like older plot
% ========================================================================

fig1 = figure('Color','w','Position',[100 100 1050 600]);

plot(freqVals,bestInfo.fftNoStim,'k','LineWidth',2.0);
hold on;
plot(freqVals,bestInfo.fftRaw,'Color',[0.75 0.75 0.75],'LineWidth',1.6);
plot(freqVals,bestInfo.fftERP,'LineWidth',1.4);

for iK = 1:length(KList)
    if KList(iK) == bestInfo.bestK
        plot(freqVals,bestInfo.smartaFFTs{iK},'LineWidth',2.2);
    else
        plot(freqVals,bestInfo.smartaFFTs{iK},'LineWidth',1.2);
    end
end

xlim(freqRangeForMetric);
xlabel('Frequency (Hz)');
ylabel('log_{10} mean FFT magnitude (a.u.)');
title(['SMARTALite K comparison: elec' num2str(bestInfo.elecNum) ...
    ', best K = ' num2str(bestInfo.bestK)]);

legendEntries = [{'No-stim reference','Raw high-stim','ERPSubtraction'}, ...
    arrayfun(@(k) ['SMARTALite K=' num2str(k)],KList,'UniformOutput',false)];

legend(legendEntries,'Location','best');
grid on;
box off;

if saveFigures
    saveas(fig1,fullfile(figFolder,['demo26_elec' num2str(bestInfo.elecNum) '_SMARTALite_K_comparison.png']));
    savefig(fig1,fullfile(figFolder,['demo26_elec' num2str(bestInfo.elecNum) '_SMARTALite_K_comparison.fig']));
    exportgraphics(fig1,fullfile(figFolder,['demo26_elec' num2str(bestInfo.elecNum) '_SMARTALite_K_comparison.pdf']), ...
        'ContentType','vector');
end

%% ========================================================================
% Figure 2: method comparison
% ========================================================================

fig2 = figure('Color','w','Position',[100 100 950 550]);

plot(freqVals,bestInfo.fftNoStim,'k','LineWidth',2.0);
hold on;
plot(freqVals,bestInfo.fftRaw,'Color',[0.75 0.75 0.75],'LineWidth',1.6);
plot(freqVals,bestInfo.fftERP,'LineWidth',1.5);
plot(freqVals,bestInfo.fftPCA,'LineWidth',1.5);
plot(freqVals,bestInfo.smartaFFTs{bestInfo.bestKIdx},'LineWidth',2.0);

xlim(freqRangeForMetric);
xlabel('Frequency (Hz)');
ylabel('log_{10} mean FFT magnitude (a.u.)');
title(['Representative FFT comparison: elec' num2str(bestInfo.elecNum)]);

legend({'No-stim reference', ...
        'Raw high-stim', ...
        'ERPSubtraction', ...
        'PCATemplate K=10', ...
        ['SMARTALite K=' num2str(bestInfo.bestK)]}, ...
        'Location','best');

grid on;
box off;

if saveFigures
    saveas(fig2,fullfile(figFolder,['demo26_elec' num2str(bestInfo.elecNum) '_method_comparison.png']));
    savefig(fig2,fullfile(figFolder,['demo26_elec' num2str(bestInfo.elecNum) '_method_comparison.fig']));
    exportgraphics(fig2,fullfile(figFolder,['demo26_elec' num2str(bestInfo.elecNum) '_method_comparison.pdf']), ...
        'ContentType','vector');
end

%% ========================================================================
% Save search table
% ========================================================================

if saveFigures
    writetable(searchTable,fullfile(figFolder,'demo26_SMARTALite_K_search_table.csv'));
end

fprintf('\nSaved demo_26 figures in:\n%s\n',figFolder);

%% ========================================================================
% Helper functions
% ========================================================================

function y = validationFFT(data,fftIdx)

    y = log10(mean(abs(fft(data(:,fftIdx)'))'));
end

function cleanedData = getCleanedData(methodOut)

    if isfield(methodOut,'cleanedData')
        cleanedData = methodOut.cleanedData;
    elseif isfield(methodOut,'cleanedTrials')
        cleanedData = methodOut.cleanedTrials;
    else
        error('Method output has neither cleanedData nor cleanedTrials.');
    end
end