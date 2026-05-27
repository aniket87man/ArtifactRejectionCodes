clear; close all; clc;

%% demo_18_SMARTAShrinkLiteAcrossGoodV1
%
% Goal:
% Test modified SMARTA-shrinkage method across the same 31 good V1 electrodes.
%
% Main question:
% Does SVD-shrinkage-based neighbor selection improve over SMARTALite K=3?
%
% Methods:
% 1. PCATemplate K=20
% 2. SMARTALite K=3
% 3. SMARTALite K=5
% 4. SMARTAShrinkLite K=3, hp=0
% 5. SMARTAShrinkLite K=3, hp=150
% 6. SMARTAShrinkLite K=5, hp=0
% 7. SMARTAShrinkLite K=5, hp=150
%
% Important:
% SMARTAShrinkLite uses shrinkage only for KNN neighbor selection.
% The final template is still the median of original stimulation-cycle segments.

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

maxElectrodesToRun = Inf;

% FFT metric settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% PCA settings
artifactWindow = [-0.02 0.4];

% SMARTA full-cycle settings from old successful version
pulseTimes = 0 + (0:7)*0.05;   % 0, 0.05, ..., 0.35 s
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Shrinkage settings
maxRank = 10;
shrinkStrength = 1.0;
useSoftShrinkage = false;

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

fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));
fprintf('FFT window samples = %d\n',length(fftIdx));
fprintf('FFT window duration = %.4f s\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% ========================================================================
% Check required functions
% ========================================================================

requiredMethods = { ...
    'PCATemplate', ...
    'SMARTALite', ...
    'SMARTAShrinkLite', ...
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
else
    warning('Could not find impedanceValues.mat. Using no impedance-based bad channels.');
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

fprintf('\nGood V1 electrodes before stim-electrode exclusion:\n');
disp(goodV1Electrodes);
fprintf('Number of good V1 electrodes = %d\n',length(goodV1Electrodes));

if excludeStimElectrode
    goodV1Electrodes = setdiff(goodV1Electrodes,stimElectrode);
    fprintf('\nExcluding stimulation electrode elec%d.\n',stimElectrode);
end

fprintf('\nFinal good V1 electrodes used in demo_18:\n');
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

% PCATemplate K=20
pcaParamsK20 = struct();
pcaParamsK20.artifactWindow = artifactWindow;
pcaParamsK20.numComponents = 20;
pcaParamsK20.removeMeanTemplate = true;
pcaParamsK20.taperEdgeMS = 2;
pcaParamsK20.doBaselineCorrection = false;

% Method list for reporting
methodNames = { ...
    'PCATemplate_K20', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTAShrink_K3_HP0', ...
    'SMARTAShrink_K3_HP150', ...
    'SMARTAShrink_K5_HP0', ...
    'SMARTAShrink_K5_HP150'};

%% ========================================================================
% Output folders
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo18_SMARTAShrinkLiteGoodV1');

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

        %% PCATemplate K=20
        pcaK20Out = PCATemplate(dataHighStim,timeVals,pcaParamsK20);
        fftPCAK20 = compute_fft_summary(pcaK20Out.cleanedData,timeVals,fftIdx);

        pcaK20Error = norm(fftPCAK20.logMeanMagnitude(freqMask) - ...
                           fftNoStim.logMeanMagnitude(freqMask));

        pcaK20Improvement = 100 * (rawError - pcaK20Error) / rawError;

        %% SMARTALite K=3
        smartaK3Opts = struct();
        smartaK3Opts.pulseTimes = pulseTimes;
        smartaK3Opts.stiFreq = smartaStimFreq;
        smartaK3Opts.prePulse = smartaPrePulse;
        smartaK3Opts.K = 3;
        smartaK3Opts.window = fftWindow;
        smartaK3Opts.computeFFT = false;

        smartaK3Out = SMARTALite(dataHighStim,timeVals,smartaK3Opts);
        smartaK3Cleaned = getCleanedData(smartaK3Out);
        fftSMARTAK3 = compute_fft_summary(smartaK3Cleaned,timeVals,fftIdx);

        smartaK3Error = norm(fftSMARTAK3.logMeanMagnitude(freqMask) - ...
                             fftNoStim.logMeanMagnitude(freqMask));

        smartaK3Improvement = 100 * (rawError - smartaK3Error) / rawError;

        %% SMARTALite K=5
        smartaK5Opts = smartaK3Opts;
        smartaK5Opts.K = 5;

        smartaK5Out = SMARTALite(dataHighStim,timeVals,smartaK5Opts);
        smartaK5Cleaned = getCleanedData(smartaK5Out);
        fftSMARTAK5 = compute_fft_summary(smartaK5Cleaned,timeVals,fftIdx);

        smartaK5Error = norm(fftSMARTAK5.logMeanMagnitude(freqMask) - ...
                             fftNoStim.logMeanMagnitude(freqMask));

        smartaK5Improvement = 100 * (rawError - smartaK5Error) / rawError;

        %% SMARTAShrinkLite variants

        shrinkErrors = nan(1,4);
        shrinkImprovements = nan(1,4);
        shrinkRanks = nan(1,4);

        shrinkConfigs = [ ...
            3   0; ...
            3 150; ...
            5   0; ...
            5 150];

        for iCfg = 1:size(shrinkConfigs,1)

            Kval = shrinkConfigs(iCfg,1);
            hpVal = shrinkConfigs(iCfg,2);

            shrinkOpts = struct();
            shrinkOpts.pulseTimes = pulseTimes;
            shrinkOpts.stiFreq = smartaStimFreq;
            shrinkOpts.prePulse = smartaPrePulse;
            shrinkOpts.K = Kval;
            shrinkOpts.hpCutoff = hpVal;
            shrinkOpts.maxRank = maxRank;
            shrinkOpts.shrinkStrength = shrinkStrength;
            shrinkOpts.useSoftShrinkage = useSoftShrinkage;
            shrinkOpts.computeFFT = false;
            shrinkOpts.window = fftWindow;

            shrinkOut = SMARTAShrinkLite(dataHighStim,timeVals,shrinkOpts);
            shrinkCleaned = getCleanedData(shrinkOut);

            fftShrink = compute_fft_summary(shrinkCleaned,timeVals,fftIdx);

            shrinkErrors(iCfg) = norm(fftShrink.logMeanMagnitude(freqMask) - ...
                                      fftNoStim.logMeanMagnitude(freqMask));

            shrinkImprovements(iCfg) = 100 * (rawError - shrinkErrors(iCfg)) / rawError;

            if isfield(shrinkOut,'diagnostics') && isfield(shrinkOut.diagnostics,'rankUse')
                shrinkRanks(iCfg) = shrinkOut.diagnostics.rankUse;
            end
        end

        %% Store row
        rowCounter = rowCounter + 1;

        resultRows(rowCounter).electrode = elecNum;
        resultRows(rowCounter).rawError = rawError;

        resultRows(rowCounter).pcaK20Error = pcaK20Error;
        resultRows(rowCounter).pcaK20Improvement = pcaK20Improvement;

        resultRows(rowCounter).smartaK3Error = smartaK3Error;
        resultRows(rowCounter).smartaK3Improvement = smartaK3Improvement;

        resultRows(rowCounter).smartaK5Error = smartaK5Error;
        resultRows(rowCounter).smartaK5Improvement = smartaK5Improvement;

        resultRows(rowCounter).shrinkK3HP0Error = shrinkErrors(1);
        resultRows(rowCounter).shrinkK3HP0Improvement = shrinkImprovements(1);
        resultRows(rowCounter).shrinkK3HP0Rank = shrinkRanks(1);

        resultRows(rowCounter).shrinkK3HP150Error = shrinkErrors(2);
        resultRows(rowCounter).shrinkK3HP150Improvement = shrinkImprovements(2);
        resultRows(rowCounter).shrinkK3HP150Rank = shrinkRanks(2);

        resultRows(rowCounter).shrinkK5HP0Error = shrinkErrors(3);
        resultRows(rowCounter).shrinkK5HP0Improvement = shrinkImprovements(3);
        resultRows(rowCounter).shrinkK5HP0Rank = shrinkRanks(3);

        resultRows(rowCounter).shrinkK5HP150Error = shrinkErrors(4);
        resultRows(rowCounter).shrinkK5HP150Improvement = shrinkImprovements(4);
        resultRows(rowCounter).shrinkK5HP150Rank = shrinkRanks(4);

        %% Print
        fprintf('Raw error                  = %.4f\n',rawError);
        fprintf('PCATemplate K=20           = %.4f, improvement %.2f %%\n',pcaK20Error,pcaK20Improvement);
        fprintf('SMARTALite K=3             = %.4f, improvement %.2f %%\n',smartaK3Error,smartaK3Improvement);
        fprintf('SMARTALite K=5             = %.4f, improvement %.2f %%\n',smartaK5Error,smartaK5Improvement);
        fprintf('SMARTAShrink K=3 HP=0      = %.4f, improvement %.2f %%, rank %d\n', ...
            shrinkErrors(1),shrinkImprovements(1),round(shrinkRanks(1)));
        fprintf('SMARTAShrink K=3 HP=150    = %.4f, improvement %.2f %%, rank %d\n', ...
            shrinkErrors(2),shrinkImprovements(2),round(shrinkRanks(2)));
        fprintf('SMARTAShrink K=5 HP=0      = %.4f, improvement %.2f %%, rank %d\n', ...
            shrinkErrors(3),shrinkImprovements(3),round(shrinkRanks(3)));
        fprintf('SMARTAShrink K=5 HP=150    = %.4f, improvement %.2f %%, rank %d\n', ...
            shrinkErrors(4),shrinkImprovements(4),round(shrinkRanks(4)));

    catch ME

        warning('Failed on elec%d: %s',elecNum,ME.message);

    end
end

if isempty(resultRows)
    error('No electrodes were processed successfully.');
end

resultsTable = struct2table(resultRows);

fprintf('\nDemo 18 SMARTAShrinkLite results table:\n');
disp(resultsTable);

%% ========================================================================
% Build matrices
% ========================================================================

errorMatrix = [ ...
    resultsTable.pcaK20Error, ...
    resultsTable.smartaK3Error, ...
    resultsTable.smartaK5Error, ...
    resultsTable.shrinkK3HP0Error, ...
    resultsTable.shrinkK3HP150Error, ...
    resultsTable.shrinkK5HP0Error, ...
    resultsTable.shrinkK5HP150Error];

improvementMatrix = [ ...
    resultsTable.pcaK20Improvement, ...
    resultsTable.smartaK3Improvement, ...
    resultsTable.smartaK5Improvement, ...
    resultsTable.shrinkK3HP0Improvement, ...
    resultsTable.shrinkK3HP150Improvement, ...
    resultsTable.shrinkK5HP0Improvement, ...
    resultsTable.shrinkK5HP150Improvement];

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
    'VariableNames',{'method','meanError','stdError','meanImprovement','stdImprovement'} ...
);

summaryTable = sortrows(summaryTable,'meanError','ascend');

fprintf('\nDemo 18 SMARTAShrinkLite summary table:\n');
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
% Figures
% ========================================================================

figure('Name','demo18 mean improvement');

bar(meanImprovement);
hold on;
errorbar(1:nMethods,meanImprovement,stdImprovement,'k.','LineWidth',1.2);

xticks(1:nMethods);
xticklabels(methodNames);
xtickangle(30);

ylabel('Improvement relative to raw high-stim (%)');
title('demo18: SMARTAShrinkLite comparison across good V1 electrodes');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo18_mean_improvement.png'));
    savefig(gcf,fullfile(figFolder,'demo18_mean_improvement.fig'));
end

figure('Name','demo18 error distribution');

boxplot(errorMatrix,'Labels',methodNames);
ylabel('FFT error vs no-stim');
title('demo18: Error distribution across good V1 electrodes');
xtickangle(30);
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo18_error_distribution.png'));
    savefig(gcf,fullfile(figFolder,'demo18_error_distribution.fig'));
end

figure('Name','demo18 electrode-wise best comparison');

plot(resultsTable.electrode,resultsTable.smartaK3Improvement,'o-','LineWidth',1.2);
hold on;
plot(resultsTable.electrode,resultsTable.shrinkK3HP0Improvement,'o-','LineWidth',1.2);
plot(resultsTable.electrode,resultsTable.shrinkK3HP150Improvement,'o-','LineWidth',1.2);
plot(resultsTable.electrode,resultsTable.pcaK20Improvement,'o-','LineWidth',1.2);

xlabel('Electrode');
ylabel('Improvement relative to raw high-stim (%)');
title('SMARTALite vs SMARTAShrinkLite vs PCATemplate');
legend({'SMARTALite K=3','Shrink K=3 HP=0','Shrink K=3 HP=150','PCATemplate K=20'}, ...
    'Location','best');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo18_electrodewise_comparison.png'));
    savefig(gcf,fullfile(figFolder,'demo18_electrodewise_comparison.fig'));
end

figure('Name','demo18 scatter SMARTALite vs best shrink');

% Pick best shrink method by mean error
shrinkMeanErrors = [ ...
    mean(resultsTable.shrinkK3HP0Error), ...
    mean(resultsTable.shrinkK3HP150Error), ...
    mean(resultsTable.shrinkK5HP0Error), ...
    mean(resultsTable.shrinkK5HP150Error)];

[~,bestShrinkIdx] = min(shrinkMeanErrors);

switch bestShrinkIdx
    case 1
        bestShrinkImp = resultsTable.shrinkK3HP0Improvement;
        bestShrinkLabel = 'Shrink K=3 HP=0';
    case 2
        bestShrinkImp = resultsTable.shrinkK3HP150Improvement;
        bestShrinkLabel = 'Shrink K=3 HP=150';
    case 3
        bestShrinkImp = resultsTable.shrinkK5HP0Improvement;
        bestShrinkLabel = 'Shrink K=5 HP=0';
    case 4
        bestShrinkImp = resultsTable.shrinkK5HP150Improvement;
        bestShrinkLabel = 'Shrink K=5 HP=150';
end

scatter(resultsTable.smartaK3Improvement,bestShrinkImp,60,resultsTable.electrode,'filled');
colorbar;

xlabel('SMARTALite K=3 improvement (%)');
ylabel([bestShrinkLabel ' improvement (%)']);
title('SMARTALite K=3 vs best SMARTAShrinkLite variant');
grid on;

hold on;
minVal = min([resultsTable.smartaK3Improvement; bestShrinkImp]);
maxVal = max([resultsTable.smartaK3Improvement; bestShrinkImp]);
plot([minVal maxVal],[minVal maxVal],'k--','LineWidth',1.2);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo18_scatter_smarta_vs_best_shrink.png'));
    savefig(gcf,fullfile(figFolder,'demo18_scatter_smarta_vs_best_shrink.fig'));
end

%% ========================================================================
% Save results
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo18_SMARTAShrinkLite_goodV1_results_' timestampString '.mat']);
resultsCsvFile = fullfile(resultsFolder,['demo18_SMARTAShrinkLite_goodV1_results_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo18_SMARTAShrinkLite_goodV1_summary_' timestampString '.csv']);

save(resultsMatFile, ...
    'resultsTable', ...
    'summaryTable', ...
    'bestMethodTable', ...
    'methodNames', ...
    'goodV1Electrodes', ...
    'pulseTimes', ...
    'smartaStimFreq', ...
    'smartaPrePulse', ...
    'maxRank', ...
    'shrinkStrength', ...
    'useSoftShrinkage', ...
    'artifactWindow', ...
    'fftWindow', ...
    'freqRangeForMetric');

writetable(resultsTable,resultsCsvFile);
writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved demo_18 results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',resultsCsvFile);
fprintf('%s\n',summaryCsvFile);

if saveFigures
    fprintf('\nSaved demo_18 figures in:\n');
    fprintf('%s\n',figFolder);
end

%% ========================================================================
% Final interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_18 SMARTAShrinkLite across good V1 complete\n');
fprintf('============================================================\n');

fprintf('Number of electrodes processed = %d\n',height(resultsTable));

fprintf('\nMethod ranking by mean FFT error:\n');
disp(summaryTable);

fprintf('\nMain comparison:\n');

meanSMARTALiteK3 = mean(resultsTable.smartaK3Error);
meanShrinkK3HP0 = mean(resultsTable.shrinkK3HP0Error);
meanShrinkK3HP150 = mean(resultsTable.shrinkK3HP150Error);

fprintf('Mean error SMARTALite K=3          = %.4f\n',meanSMARTALiteK3);
fprintf('Mean error SMARTAShrink K=3 HP=0   = %.4f\n',meanShrinkK3HP0);
fprintf('Mean error SMARTAShrink K=3 HP=150 = %.4f\n',meanShrinkK3HP150);

if meanShrinkK3HP0 < meanSMARTALiteK3 || meanShrinkK3HP150 < meanSMARTALiteK3
    fprintf('At least one SMARTAShrinkLite variant improved over SMARTALite K=3.\n');
else
    fprintf('SMARTALite K=3 remains better than the tested SMARTAShrinkLite variants.\n');
end

%% ========================================================================
% Local helper
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