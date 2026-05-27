clear; close all; clc;

%% demo_17_SMARTALiteAcrossGoodElectrodes
%
% Goal:
% Test SMARTALite across the same corrected good V1 electrode set used in demo_14.
%
% Main question:
% Does SMARTALite K=5 outperform PCATemplate K=20 across all good V1 electrodes?
%
% Methods:
% 1. ERPSubtraction
% 2. PCATemplate K=3
% 3. PCATemplate K=20
% 4. SMARTALite K=3
% 5. SMARTALite K=5
% 6. SMARTALite K=10
% 7. SMARTALite K=20
%
% Important:
% SMARTALite K is the number of nearest-neighbor stimulation-cycle segments.
% It is NOT the same as PCATemplate K, which is number of PCA components removed.

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

% Electrode selection
v1Electrodes = 1:48;
stimElectrode = 1;
excludeStimElectrode = true;

% Run settings
maxElectrodesToRun = Inf;   % Use Inf for all 31 good V1 electrodes

% FFT metric settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% PCA settings
artifactWindow = [-0.02 0.4];

% SMARTALite old successful full-cycle settings
pulseTimes = 0 + (0:7)*0.05;   % 0, 0.05, ..., 0.35 s
smartaStimFreq = 20;           % 50 ms stimulation cycle
smartaPrePulse = 0.0005;       % 0.5 ms before pulse

% Test several neighbor counts
smartaKList = [3 5 10 20];

% Save results and figures
saveFigures = true;

%% ========================================================================
% Paths and loading
% ========================================================================

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder   = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

if ~exist(lfpInfoFile,'file')
    error('Could not find lfpInfo file: %s',lfpInfoFile);
end

if ~exist(paramFile,'file')
    error('Could not find parameterCombinations file: %s',paramFile);
end

I = load(lfpInfoFile);
P = load(paramFile);

timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

noStimTrials = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));
fprintf('FFT window samples = %d\n',length(fftIdx));
fprintf('FFT window duration = %.4f s\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% ========================================================================
% Check required functions
% ========================================================================

requiredMethods = { ...
    'ERPSubtraction', ...
    'PCATemplate', ...
    'SMARTALite', ...
    'compute_fft_summary'};

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
% Corrected good V1 electrode selection
% ========================================================================

badChannels = [];

impedanceFileName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat');

if exist(impedanceFileName,'file')
    Z = load(impedanceFileName);

    if isfield(Z,'impedanceValues')
        impedanceValues = Z.impedanceValues;
    elseif isfield(Z,'electrodeImpedances')
        impedanceValues = Z.electrodeImpedances;
    else
        error('Could not find impedanceValues or electrodeImpedances in impedance file.');
    end

    badImpedanceCutoff = 2500;
    badChannels = unique([find(impedanceValues > badImpedanceCutoff), find(isnan(impedanceValues))]);

    fprintf('\nImpedance file used:\n');
    disp(impedanceFileName);
else
    warning('Could not find impedanceValues.mat. Using no impedance-based bad channels.');
    badChannels = [];
end

rfDataFileName = [subjectName gridType 'RFData.mat'];
rfDataPath = which(rfDataFileName);

if isempty(rfDataPath)
    error('Could not find %s on MATLAB path.',rfDataFileName);
end

rfData = load(rfDataPath);

if ~isfield(rfData,'highRMSElectrodes')
    error('RFData file does not contain highRMSElectrodes.');
end

highRMSElectrodes = rfData.highRMSElectrodes(:)';

goodElectrodesAll = setdiff(highRMSElectrodes,badChannels);
goodV1Electrodes = intersect(goodElectrodesAll,v1Electrodes);

fprintf('\nRFData file used:\n');
disp(rfDataPath);

fprintf('\nGood V1 electrodes before stim-electrode exclusion:\n');
disp(goodV1Electrodes);
fprintf('Number of good V1 electrodes = %d\n',length(goodV1Electrodes));

if excludeStimElectrode
    goodV1Electrodes = setdiff(goodV1Electrodes,stimElectrode);
    fprintf('\nExcluding stimulation electrode elec%d.\n',stimElectrode);
end

fprintf('\nFinal good V1 electrodes used in demo_17:\n');
disp(goodV1Electrodes);
fprintf('Final number of electrodes = %d\n',length(goodV1Electrodes));

%% ========================================================================
% Find electrode files
% ========================================================================

elecFiles = dir(fullfile(lfpFolder,'elec*.mat'));

if isempty(elecFiles)
    error('No elec*.mat files found in %s',lfpFolder);
end

elecNums = nan(length(elecFiles),1);

for i = 1:length(elecFiles)
    token = regexp(elecFiles(i).name,'^elec(\d+)\.mat$','tokens');

    if ~isempty(token)
        elecNums(i) = str2double(token{1}{1});
    end
end

validFileIdx = ~isnan(elecNums);
elecFiles = elecFiles(validFileIdx);
elecNums = elecNums(validFileIdx);

[elecNums,sortIdx] = sort(elecNums);
elecFiles = elecFiles(sortIdx);

keepIdx = ismember(elecNums,goodV1Electrodes);

elecNums = elecNums(keepIdx);
elecFiles = elecFiles(keepIdx);

if isempty(elecNums)
    error('No electrode files matched the corrected good V1 electrode list.');
end

if isfinite(maxElectrodesToRun)
    nToRun = min(maxElectrodesToRun,length(elecNums));
    elecNums = elecNums(1:nToRun);
    elecFiles = elecFiles(1:nToRun);
end

fprintf('\nNumber of electrode files selected for analysis: %d\n',length(elecNums));
fprintf('Electrode files selected:\n');
disp(elecNums');

%% ========================================================================
% Method parameters
% ========================================================================

% ERPSubtraction
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

% PCATemplate K=3
pcaParamsK3 = struct();
pcaParamsK3.artifactWindow = artifactWindow;
pcaParamsK3.numComponents = 3;
pcaParamsK3.removeMeanTemplate = true;
pcaParamsK3.taperEdgeMS = 2;
pcaParamsK3.doBaselineCorrection = false;

% PCATemplate K=20
pcaParamsK20 = struct();
pcaParamsK20.artifactWindow = artifactWindow;
pcaParamsK20.numComponents = 20;
pcaParamsK20.removeMeanTemplate = true;
pcaParamsK20.taperEdgeMS = 2;
pcaParamsK20.doBaselineCorrection = false;

%% ========================================================================
% Output folders
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo17_SMARTALiteGoodV1');

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Loop across electrodes
% ========================================================================

resultRows = struct([]);
rowCounter = 0;

for iElec = 1:length(elecFiles)

    elecNum = elecNums(iElec);
    elecFile = fullfile(lfpFolder,elecFiles(iElec).name);

    fprintf('\n============================================================\n');
    fprintf('Processing electrode %d/%d: elec%d\n',iElec,length(elecFiles),elecNum);
    fprintf('============================================================\n');

    try
        D = load(elecFile);

        if ~isfield(D,'analogData')
            warning('Skipping elec%d: analogData not found.',elecNum);
            continue;
        end

        analogData = D.analogData;

        dataNoStim = analogData(noStimTrials,:);
        dataHighStim = analogData(highStimTrials,:);

        %% Reference FFT
        fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);
        fftHighRaw = compute_fft_summary(dataHighStim,timeVals,fftIdx);

        freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
                   fftNoStim.freqAxis <= freqRangeForMetric(2);

        rawError = norm(fftHighRaw.logMeanMagnitude(freqMask) - ...
                        fftNoStim.logMeanMagnitude(freqMask));

        %% ERPSubtraction
        erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);
        fftERP = compute_fft_summary(erpOut.cleanedData,timeVals,fftIdx);

        erpError = norm(fftERP.logMeanMagnitude(freqMask) - ...
                        fftNoStim.logMeanMagnitude(freqMask));

        erpImprovement = 100 * (rawError - erpError) / rawError;

        %% PCATemplate K=3
        pcaK3Out = PCATemplate(dataHighStim,timeVals,pcaParamsK3);
        fftPCAK3 = compute_fft_summary(pcaK3Out.cleanedData,timeVals,fftIdx);

        pcaK3Error = norm(fftPCAK3.logMeanMagnitude(freqMask) - ...
                          fftNoStim.logMeanMagnitude(freqMask));

        pcaK3Improvement = 100 * (rawError - pcaK3Error) / rawError;

        %% PCATemplate K=20
        pcaK20Out = PCATemplate(dataHighStim,timeVals,pcaParamsK20);
        fftPCAK20 = compute_fft_summary(pcaK20Out.cleanedData,timeVals,fftIdx);

        pcaK20Error = norm(fftPCAK20.logMeanMagnitude(freqMask) - ...
                           fftNoStim.logMeanMagnitude(freqMask));

        pcaK20Improvement = 100 * (rawError - pcaK20Error) / rawError;

        %% SMARTALite K sweep
        smartaErrors = nan(1,length(smartaKList));
        smartaImprovements = nan(1,length(smartaKList));

        for iK = 1:length(smartaKList)

            K = smartaKList(iK);

            smartaOpts = struct();
            smartaOpts.pulseTimes = pulseTimes;
            smartaOpts.stiFreq = smartaStimFreq;
            smartaOpts.prePulse = smartaPrePulse;
            smartaOpts.K = K;
            smartaOpts.window = fftWindow;
            smartaOpts.computeFFT = false;

            try
                smartaOut = SMARTALite(dataHighStim,timeVals,smartaOpts);

                if isfield(smartaOut,'cleanedData')
                    smartaCleaned = smartaOut.cleanedData;
                elseif isfield(smartaOut,'cleanedTrials')
                    smartaCleaned = smartaOut.cleanedTrials;
                else
                    error('SMARTALite output has neither cleanedData nor cleanedTrials.');
                end

                fftSMARTA = compute_fft_summary(smartaCleaned,timeVals,fftIdx);

                smartaErrors(iK) = norm(fftSMARTA.logMeanMagnitude(freqMask) - ...
                                        fftNoStim.logMeanMagnitude(freqMask));

                smartaImprovements(iK) = 100 * (rawError - smartaErrors(iK)) / rawError;

            catch ME
                warning('SMARTALite K=%d failed on elec%d: %s',K,elecNum,ME.message);
            end
        end

        %% Store row
        rowCounter = rowCounter + 1;

        resultRows(rowCounter).electrode = elecNum;
        resultRows(rowCounter).rawError = rawError;

        resultRows(rowCounter).erpError = erpError;
        resultRows(rowCounter).erpImprovement = erpImprovement;

        resultRows(rowCounter).pcaK3Error = pcaK3Error;
        resultRows(rowCounter).pcaK3Improvement = pcaK3Improvement;

        resultRows(rowCounter).pcaK20Error = pcaK20Error;
        resultRows(rowCounter).pcaK20Improvement = pcaK20Improvement;

        for iK = 1:length(smartaKList)
            K = smartaKList(iK);

            resultRows(rowCounter).(['smartaK' num2str(K) 'Error']) = smartaErrors(iK);
            resultRows(rowCounter).(['smartaK' num2str(K) 'Improvement']) = smartaImprovements(iK);
        end

        %% Print
        fprintf('Raw error          = %.4f\n',rawError);
        fprintf('ERPSubtraction     = %.4f, improvement %.2f %%\n',erpError,erpImprovement);
        fprintf('PCATemplate K=3    = %.4f, improvement %.2f %%\n',pcaK3Error,pcaK3Improvement);
        fprintf('PCATemplate K=20   = %.4f, improvement %.2f %%\n',pcaK20Error,pcaK20Improvement);

        for iK = 1:length(smartaKList)
            K = smartaKList(iK);
            fprintf('SMARTALite K=%-3d  = %.4f, improvement %.2f %%\n', ...
                K,smartaErrors(iK),smartaImprovements(iK));
        end

    catch ME

        warning('Failed on elec%d: %s',elecNum,ME.message);

    end
end

if isempty(resultRows)
    error('No electrodes were processed successfully.');
end

resultsTable = struct2table(resultRows);

fprintf('\nDemo 17 SMARTALite across-good-V1 results table:\n');
disp(resultsTable);

%% ========================================================================
% Build method matrices
% ========================================================================

methodNames = {'ERPSubtraction','PCATemplate_K3','PCATemplate_K20'};

errorMatrix = [ ...
    resultsTable.erpError, ...
    resultsTable.pcaK3Error, ...
    resultsTable.pcaK20Error];

improvementMatrix = [ ...
    resultsTable.erpImprovement, ...
    resultsTable.pcaK3Improvement, ...
    resultsTable.pcaK20Improvement];

for iK = 1:length(smartaKList)
    K = smartaKList(iK);

    methodNames{end+1} = ['SMARTALite_K' num2str(K)];

    errorMatrix = [errorMatrix, resultsTable.(['smartaK' num2str(K) 'Error'])];
    improvementMatrix = [improvementMatrix, resultsTable.(['smartaK' num2str(K) 'Improvement'])];
end

%% ========================================================================
% Summary statistics
% ========================================================================

nMethods = length(methodNames);

meanError = nan(1,nMethods);
stdError = nan(1,nMethods);
meanImprovement = nan(1,nMethods);
stdImprovement = nan(1,nMethods);

for iMethod = 1:nMethods
    xErr = errorMatrix(:,iMethod);
    xImp = improvementMatrix(:,iMethod);

    xErr = xErr(~isnan(xErr));
    xImp = xImp(~isnan(xImp));

    meanError(iMethod) = mean(xErr);
    stdError(iMethod) = std(xErr);

    meanImprovement(iMethod) = mean(xImp);
    stdImprovement(iMethod) = std(xImp);
end

summaryTable = table( ...
    methodNames(:), ...
    meanError(:), ...
    stdError(:), ...
    meanImprovement(:), ...
    stdImprovement(:), ...
    'VariableNames', { ...
        'method', ...
        'meanError', ...
        'stdError', ...
        'meanImprovement', ...
        'stdImprovement' ...
    });

summaryTable = sortrows(summaryTable,'meanError','ascend');

fprintf('\nDemo 17 SMARTALite across-good-V1 summary table:\n');
disp(summaryTable);

%% ========================================================================
% Best method per electrode
% ========================================================================

errorMatrixForBest = errorMatrix;
errorMatrixForBest(isnan(errorMatrixForBest)) = inf;

[~,bestMethodIdx] = min(errorMatrixForBest,[],2);

bestMethodPerElectrode = methodNames(bestMethodIdx)';

bestMethodTable = table( ...
    resultsTable.electrode, ...
    bestMethodPerElectrode, ...
    'VariableNames',{'electrode','bestMethod'} ...
);

fprintf('\nBest method per electrode:\n');
disp(bestMethodTable);

fprintf('\nBest method counts:\n');

for iMethod = 1:length(methodNames)
    nBest = sum(bestMethodIdx == iMethod);
    fprintf('%s: %d electrodes\n',methodNames{iMethod},nBest);
end

%% ========================================================================
% Figure 1: Mean improvement
% ========================================================================

figure('Name','demo17 SMARTALite across good V1: mean improvement');

bar(meanImprovement);
hold on;
errorbar(1:nMethods,meanImprovement,stdImprovement,'k.','LineWidth',1.2);

xticks(1:nMethods);
xticklabels(methodNames);
xtickangle(30);

ylabel('Improvement relative to raw high-stim (%)');
title('demo17: SMARTALite across good V1 electrodes');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo17_mean_improvement.png'));
    savefig(gcf,fullfile(figFolder,'demo17_mean_improvement.fig'));
end

%% ========================================================================
% Figure 2: Error distribution
% ========================================================================

figure('Name','demo17 SMARTALite across good V1: error distribution');

boxplot(errorMatrix,'Labels',methodNames);
ylabel('FFT error vs no-stim');
title('demo17: Error distribution across good V1 electrodes');
xtickangle(30);
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo17_error_distribution.png'));
    savefig(gcf,fullfile(figFolder,'demo17_error_distribution.fig'));
end

%% ========================================================================
% Figure 3: Electrode-wise comparison for main methods
% ========================================================================

figure('Name','demo17 SMARTALite across good V1: electrode-wise main comparison');

plot(resultsTable.electrode,resultsTable.pcaK20Improvement,'o-','LineWidth',1.2);
hold on;
plot(resultsTable.electrode,resultsTable.smartaK5Improvement,'o-','LineWidth',1.2);

xlabel('Electrode');
ylabel('Improvement relative to raw high-stim (%)');
title('PCATemplate K=20 vs SMARTALite K=5');
legend({'PCATemplate K=20','SMARTALite K=5'},'Location','best');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo17_pcaK20_vs_smartaK5_electrodewise.png'));
    savefig(gcf,fullfile(figFolder,'demo17_pcaK20_vs_smartaK5_electrodewise.fig'));
end

%% ========================================================================
% Figure 4: SMARTALite K sweep
% ========================================================================

smartaMeanError = nan(1,length(smartaKList));
smartaStdError = nan(1,length(smartaKList));

for iK = 1:length(smartaKList)
    K = smartaKList(iK);
    x = resultsTable.(['smartaK' num2str(K) 'Error']);
    x = x(~isnan(x));
    smartaMeanError(iK) = mean(x);
    smartaStdError(iK) = std(x);
end

figure('Name','demo17 SMARTALite K sweep');

errorbar(smartaKList,smartaMeanError,smartaStdError,'o-','LineWidth',1.2);
xlabel('SMARTALite K');
ylabel('Mean FFT error vs no-stim');
title('SMARTALite neighbor-count sweep across good V1 electrodes');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo17_SMARTALite_K_sweep.png'));
    savefig(gcf,fullfile(figFolder,'demo17_SMARTALite_K_sweep.fig'));
end

%% ========================================================================
% Figure 5: Scatter PCATemplate K=20 vs SMARTALite K=5
% ========================================================================

figure('Name','demo17 PCATemplate K20 vs SMARTALite K5');

scatter(resultsTable.pcaK20Improvement,resultsTable.smartaK5Improvement,60,resultsTable.electrode,'filled');
colorbar;

xlabel('PCATemplate K=20 improvement (%)');
ylabel('SMARTALite K=5 improvement (%)');
title('SMARTALite K=5 vs PCATemplate K=20 across good V1 electrodes');
grid on;

hold on;
minVal = min([resultsTable.pcaK20Improvement; resultsTable.smartaK5Improvement]);
maxVal = max([resultsTable.pcaK20Improvement; resultsTable.smartaK5Improvement]);
plot([minVal maxVal],[minVal maxVal],'k--','LineWidth',1.2);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo17_scatter_pcaK20_vs_smartaK5.png'));
    savefig(gcf,fullfile(figFolder,'demo17_scatter_pcaK20_vs_smartaK5.fig'));
end

%% ========================================================================
% Save results
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo17_SMARTALite_goodV1_results_' timestampString '.mat']);
resultsCsvFile = fullfile(resultsFolder,['demo17_SMARTALite_goodV1_results_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo17_SMARTALite_goodV1_summary_' timestampString '.csv']);

save(resultsMatFile, ...
    'resultsTable', ...
    'summaryTable', ...
    'bestMethodTable', ...
    'methodNames', ...
    'goodV1Electrodes', ...
    'goodElectrodesAll', ...
    'highRMSElectrodes', ...
    'badChannels', ...
    'rfDataPath', ...
    'pulseTimes', ...
    'smartaKList', ...
    'smartaStimFreq', ...
    'smartaPrePulse', ...
    'artifactWindow', ...
    'fftWindow', ...
    'freqRangeForMetric');

writetable(resultsTable,resultsCsvFile);
writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved demo_17 results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',resultsCsvFile);
fprintf('%s\n',summaryCsvFile);

if saveFigures
    fprintf('\nSaved demo_17 figures in:\n');
    fprintf('%s\n',figFolder);
end

%% ========================================================================
% Final interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_17 SMARTALite across good V1 complete\n');
fprintf('============================================================\n');

fprintf('Number of electrodes processed = %d\n',height(resultsTable));
fprintf('Pulse times used:\n');
disp(pulseTimes);

fprintf('\nMethod ranking by mean FFT error:\n');
disp(summaryTable);

fprintf('\nMain comparison:\n');

pcaK20Mean = mean(resultsTable.pcaK20Error);
smartaK5Mean = mean(resultsTable.smartaK5Error);

fprintf('Mean error PCATemplate K=20 = %.4f\n',pcaK20Mean);
fprintf('Mean error SMARTALite K=5   = %.4f\n',smartaK5Mean);

if smartaK5Mean < pcaK20Mean
    fprintf('SMARTALite K=5 outperformed PCATemplate K=20 on mean FFT error.\n');
else
    fprintf('PCATemplate K=20 outperformed SMARTALite K=5 on mean FFT error.\n');
end

fprintf('\nReminder:\n');
fprintf('PCATemplate K=20 removes 20 PCA components and may be aggressive.\n');
fprintf('SMARTALite K means nearest-neighbor template count, not PCA components.\n');