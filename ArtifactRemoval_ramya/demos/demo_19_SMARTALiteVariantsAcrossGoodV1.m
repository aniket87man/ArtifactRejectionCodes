clear; close all; clc;

%% demo_19_SMARTALiteVariantsAcrossGoodV1
%
% Goal:
% Test practical variants of the current best method, SMARTALite K=3.
%
% Current best from demo_17/demo_18:
%   SMARTALite K=3
%
% Variants tested:
% 1. PCATemplate K=20
% 2. SMARTALite K=3
% 3. SMARTALite K=5
% 4. SMARTALite K=3 leave-one-trial-out
% 5. SMARTALite ensemble K=3+5
% 6. SMARTALite K=3 + residual ERP alpha=0.25
% 7. SMARTALite K=3 + residual ERP alpha=0.50
% 8. SMARTALite K=3 + residual ERP alpha=1.00
%
% Important:
% The leave-one-trial-out version is a validation method:
% for each trial, it excludes all segments from that same trial
% when estimating templates.

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

% SMARTALite full-cycle settings
pulseTimes = 0 + (0:7)*0.05;   % 0, 0.05, ..., 0.35 s
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Residual ERP settings
residualWindow = [0 0.4];
residualTaperMS = 2;
residualAlphaList = [0.25 0.50 1.00];

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

fprintf('\nFinal good V1 electrodes used in demo_19:\n');
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

% SMARTALite K=3
smartaK3Opts = struct();
smartaK3Opts.pulseTimes = pulseTimes;
smartaK3Opts.stiFreq = smartaStimFreq;
smartaK3Opts.prePulse = smartaPrePulse;
smartaK3Opts.K = 3;
smartaK3Opts.window = fftWindow;
smartaK3Opts.computeFFT = false;

% SMARTALite K=5
smartaK5Opts = smartaK3Opts;
smartaK5Opts.K = 5;

% LOTO local variant
lotoOpts = smartaK3Opts;
lotoOpts.excludeSameTrial = true;

methodNames = { ...
    'PCATemplate_K20', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_K3_LOTO', ...
    'SMARTALite_Ensemble_K3K5', ...
    'SMARTALite_K3_ResidualERP_025', ...
    'SMARTALite_K3_ResidualERP_050', ...
    'SMARTALite_K3_ResidualERP_100'};

%% ========================================================================
% Output folders
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo19_SMARTALiteVariantsGoodV1');

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
        pcaK20Cleaned = getCleanedData(pcaK20Out);
        pcaK20Error = computeMetricError(pcaK20Cleaned,timeVals,fftIdx,fftNoStim,freqMask);
        pcaK20Improvement = 100 * (rawError - pcaK20Error) / rawError;

        %% SMARTALite K=3
        smartaK3Out = SMARTALite(dataHighStim,timeVals,smartaK3Opts);
        smartaK3Cleaned = getCleanedData(smartaK3Out);
        smartaK3Error = computeMetricError(smartaK3Cleaned,timeVals,fftIdx,fftNoStim,freqMask);
        smartaK3Improvement = 100 * (rawError - smartaK3Error) / rawError;

        %% SMARTALite K=5
        smartaK5Out = SMARTALite(dataHighStim,timeVals,smartaK5Opts);
        smartaK5Cleaned = getCleanedData(smartaK5Out);
        smartaK5Error = computeMetricError(smartaK5Cleaned,timeVals,fftIdx,fftNoStim,freqMask);
        smartaK5Improvement = 100 * (rawError - smartaK5Error) / rawError;

        %% SMARTALite K=3 leave-one-trial-out local variant
        lotoOut = SMARTALiteLocalFullCycle(dataHighStim,timeVals,lotoOpts);
        lotoCleaned = getCleanedData(lotoOut);
        lotoError = computeMetricError(lotoCleaned,timeVals,fftIdx,fftNoStim,freqMask);
        lotoImprovement = 100 * (rawError - lotoError) / rawError;

        %% Ensemble K=3 and K=5
        ensembleCleaned = 0.5 * (smartaK3Cleaned + smartaK5Cleaned);
        ensembleError = computeMetricError(ensembleCleaned,timeVals,fftIdx,fftNoStim,freqMask);
        ensembleImprovement = 100 * (rawError - ensembleError) / rawError;

        %% Residual ERP variants after SMARTALite K=3
        residualErrors = nan(1,length(residualAlphaList));
        residualImprovements = nan(1,length(residualAlphaList));

        for iAlpha = 1:length(residualAlphaList)

            alpha = residualAlphaList(iAlpha);

            residualCleaned = subtractResidualERP( ...
                smartaK3Cleaned,timeVals,residualWindow,alpha,residualTaperMS);

            residualErrors(iAlpha) = computeMetricError( ...
                residualCleaned,timeVals,fftIdx,fftNoStim,freqMask);

            residualImprovements(iAlpha) = 100 * (rawError - residualErrors(iAlpha)) / rawError;
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

        resultRows(rowCounter).lotoK3Error = lotoError;
        resultRows(rowCounter).lotoK3Improvement = lotoImprovement;

        resultRows(rowCounter).ensembleK3K5Error = ensembleError;
        resultRows(rowCounter).ensembleK3K5Improvement = ensembleImprovement;

        resultRows(rowCounter).resERP025Error = residualErrors(1);
        resultRows(rowCounter).resERP025Improvement = residualImprovements(1);

        resultRows(rowCounter).resERP050Error = residualErrors(2);
        resultRows(rowCounter).resERP050Improvement = residualImprovements(2);

        resultRows(rowCounter).resERP100Error = residualErrors(3);
        resultRows(rowCounter).resERP100Improvement = residualImprovements(3);

        %% Print
        fprintf('Raw error                          = %.4f\n',rawError);
        fprintf('PCATemplate K=20                   = %.4f, improvement %.2f %%\n',pcaK20Error,pcaK20Improvement);
        fprintf('SMARTALite K=3                     = %.4f, improvement %.2f %%\n',smartaK3Error,smartaK3Improvement);
        fprintf('SMARTALite K=5                     = %.4f, improvement %.2f %%\n',smartaK5Error,smartaK5Improvement);
        fprintf('SMARTALite K=3 LOTO                = %.4f, improvement %.2f %%\n',lotoError,lotoImprovement);
        fprintf('SMARTALite ensemble K3K5           = %.4f, improvement %.2f %%\n',ensembleError,ensembleImprovement);
        fprintf('SMARTALite K3 + residual ERP 0.25  = %.4f, improvement %.2f %%\n',residualErrors(1),residualImprovements(1));
        fprintf('SMARTALite K3 + residual ERP 0.50  = %.4f, improvement %.2f %%\n',residualErrors(2),residualImprovements(2));
        fprintf('SMARTALite K3 + residual ERP 1.00  = %.4f, improvement %.2f %%\n',residualErrors(3),residualImprovements(3));

    catch ME

        warning('Failed on elec%d: %s',elecNum,ME.message);

    end
end

if isempty(resultRows)
    error('No electrodes were processed successfully.');
end

resultsTable = struct2table(resultRows);

fprintf('\nDemo 19 SMARTALite variants results table:\n');
disp(resultsTable);

%% ========================================================================
% Build matrices
% ========================================================================

errorMatrix = [ ...
    resultsTable.pcaK20Error, ...
    resultsTable.smartaK3Error, ...
    resultsTable.smartaK5Error, ...
    resultsTable.lotoK3Error, ...
    resultsTable.ensembleK3K5Error, ...
    resultsTable.resERP025Error, ...
    resultsTable.resERP050Error, ...
    resultsTable.resERP100Error];

improvementMatrix = [ ...
    resultsTable.pcaK20Improvement, ...
    resultsTable.smartaK3Improvement, ...
    resultsTable.smartaK5Improvement, ...
    resultsTable.lotoK3Improvement, ...
    resultsTable.ensembleK3K5Improvement, ...
    resultsTable.resERP025Improvement, ...
    resultsTable.resERP050Improvement, ...
    resultsTable.resERP100Improvement];

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

fprintf('\nDemo 19 SMARTALite variants summary table:\n');
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

figure('Name','demo19 mean improvement');

bar(meanImprovement);
hold on;
errorbar(1:nMethods,meanImprovement,stdImprovement,'k.','LineWidth',1.2);

xticks(1:nMethods);
xticklabels(methodNames);
xtickangle(30);

ylabel('Improvement relative to raw high-stim (%)');
title('demo19: SMARTALite variants across good V1 electrodes');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo19_mean_improvement.png'));
    savefig(gcf,fullfile(figFolder,'demo19_mean_improvement.fig'));
end

figure('Name','demo19 error distribution');

boxplot(errorMatrix,'Labels',methodNames);
ylabel('FFT error vs no-stim');
title('demo19: Error distribution across good V1 electrodes');
xtickangle(30);
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo19_error_distribution.png'));
    savefig(gcf,fullfile(figFolder,'demo19_error_distribution.fig'));
end

figure('Name','demo19 electrode-wise comparison');

plot(resultsTable.electrode,resultsTable.smartaK3Improvement,'o-','LineWidth',1.2);
hold on;
plot(resultsTable.electrode,resultsTable.lotoK3Improvement,'o-','LineWidth',1.2);
plot(resultsTable.electrode,resultsTable.ensembleK3K5Improvement,'o-','LineWidth',1.2);
plot(resultsTable.electrode,resultsTable.resERP025Improvement,'o-','LineWidth',1.2);
plot(resultsTable.electrode,resultsTable.pcaK20Improvement,'o-','LineWidth',1.2);

xlabel('Electrode');
ylabel('Improvement relative to raw high-stim (%)');
title('Main SMARTALite variants by electrode');
legend({ ...
    'SMARTALite K=3', ...
    'SMARTALite K=3 LOTO', ...
    'Ensemble K3K5', ...
    'K3 + residual ERP 0.25', ...
    'PCATemplate K=20'}, ...
    'Location','best');
grid on;

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo19_electrodewise_comparison.png'));
    savefig(gcf,fullfile(figFolder,'demo19_electrodewise_comparison.fig'));
end

figure('Name','demo19 scatter baseline vs LOTO');

scatter(resultsTable.smartaK3Improvement,resultsTable.lotoK3Improvement,60,resultsTable.electrode,'filled');
colorbar;

xlabel('SMARTALite K=3 improvement (%)');
ylabel('SMARTALite K=3 LOTO improvement (%)');
title('Validation: standard SMARTALite vs leave-one-trial-out');
grid on;

hold on;
minVal = min([resultsTable.smartaK3Improvement; resultsTable.lotoK3Improvement]);
maxVal = max([resultsTable.smartaK3Improvement; resultsTable.lotoK3Improvement]);
plot([minVal maxVal],[minVal maxVal],'k--','LineWidth',1.2);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo19_scatter_smartaK3_vs_loto.png'));
    savefig(gcf,fullfile(figFolder,'demo19_scatter_smartaK3_vs_loto.fig'));
end

%% ========================================================================
% Save results
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo19_SMARTALite_variants_goodV1_results_' timestampString '.mat']);
resultsCsvFile = fullfile(resultsFolder,['demo19_SMARTALite_variants_goodV1_results_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo19_SMARTALite_variants_goodV1_summary_' timestampString '.csv']);

save(resultsMatFile, ...
    'resultsTable', ...
    'summaryTable', ...
    'bestMethodTable', ...
    'methodNames', ...
    'goodV1Electrodes', ...
    'pulseTimes', ...
    'smartaStimFreq', ...
    'smartaPrePulse', ...
    'artifactWindow', ...
    'fftWindow', ...
    'freqRangeForMetric', ...
    'residualWindow', ...
    'residualAlphaList');

writetable(resultsTable,resultsCsvFile);
writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved demo_19 results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',resultsCsvFile);
fprintf('%s\n',summaryCsvFile);

if saveFigures
    fprintf('\nSaved demo_19 figures in:\n');
    fprintf('%s\n',figFolder);
end

%% ========================================================================
% Final interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_19 SMARTALite variants across good V1 complete\n');
fprintf('============================================================\n');

fprintf('Number of electrodes processed = %d\n',height(resultsTable));

fprintf('\nMethod ranking by mean FFT error:\n');
disp(summaryTable);

fprintf('\nKey comparisons:\n');
fprintf('Mean error SMARTALite K=3      = %.4f\n',mean(resultsTable.smartaK3Error));
fprintf('Mean error SMARTALite K=3 LOTO = %.4f\n',mean(resultsTable.lotoK3Error));
fprintf('Mean error ensemble K3K5       = %.4f\n',mean(resultsTable.ensembleK3K5Error));
fprintf('Mean error residual ERP 0.25   = %.4f\n',mean(resultsTable.resERP025Error));

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

function err = computeMetricError(cleanedData,timeVals,fftIdx,fftNoStim,freqMask)

    fftCleaned = compute_fft_summary(cleanedData,timeVals,fftIdx);

    err = norm(fftCleaned.logMeanMagnitude(freqMask) - ...
               fftNoStim.logMeanMagnitude(freqMask));
end

function cleanedOut = subtractResidualERP(cleanedIn,timeVals,subtractWindow,alpha,taperMS)

    cleanedOut = cleanedIn;

    idx = timeVals >= subtractWindow(1) & timeVals < subtractWindow(2);

    if ~any(idx)
        error('No samples found in residual ERP subtractWindow.');
    end

    template = median(cleanedIn,1);

    Fs = 1/median(diff(timeVals));
    edgeSamples = round((taperMS/1000)*Fs);

    segment = template(idx);
    segment = applyEdgeTaper(segment,edgeSamples);

    cleanedOut(:,idx) = cleanedOut(:,idx) - alpha * segment;
end

function y = applyEdgeTaper(x,edgeSamples)

    y = x;

    if edgeSamples <= 0
        return;
    end

    n = length(x);
    edgeSamples = min(edgeSamples,floor(n/2));

    if edgeSamples <= 0
        return;
    end

    t = linspace(0,pi/2,edgeSamples);

    leftTaper = sin(t).^2;
    rightTaper = fliplr(leftTaper);

    y(1:edgeSamples) = y(1:edgeSamples) .* leftTaper;
    y(end-edgeSamples+1:end) = y(end-edgeSamples+1:end) .* rightTaper;
end

function out = SMARTALiteLocalFullCycle(trialData,timeVals,opts)
% Local SMARTALite implementation with leave-one-trial-out support.
%
% Uses full stimulation-cycle segments:
%   pulse time - prePulse to pulse time + 1/stiFreq

    if ~isfield(opts,'pulseTimes') || isempty(opts.pulseTimes)
        error('opts.pulseTimes must be provided.');
    end

    if ~isfield(opts,'stiFreq') || isempty(opts.stiFreq)
        error('opts.stiFreq must be provided.');
    end

    if ~isfield(opts,'prePulse') || isempty(opts.prePulse)
        opts.prePulse = 0.0005;
    end

    if ~isfield(opts,'K') || isempty(opts.K)
        opts.K = 3;
    end

    if ~isfield(opts,'excludeSameTrial') || isempty(opts.excludeSameTrial)
        opts.excludeSameTrial = false;
    end

    timeVals = timeVals(:).';

    [nTrials,nSamples] = size(trialData);
    Fs = 1/median(diff(timeVals));

    pulseTimes = opts.pulseTimes(:).';
    nPulsesOriginal = length(pulseTimes);

    stPoint = -round(opts.prePulse * Fs);
    edPoint = round((1/opts.stiFreq) * Fs);

    segLen = edPoint - stPoint + 1;

    validPulseMask = true(size(pulseTimes));

    for iPulse = 1:nPulsesOriginal
        centerIdx = nearestIndex(timeVals,pulseTimes(iPulse));

        if centerIdx + stPoint < 1 || centerIdx + edPoint > nSamples
            validPulseMask(iPulse) = false;
        end
    end

    pulseTimes = pulseTimes(validPulseMask);
    nPulses = length(pulseTimes);

    if nPulses == 0
        error('No pulseTimes remain after boundary checking.');
    end

    nSegments = nTrials * nPulses;

    X = zeros(nSegments,segLen);
    segmentTrial = zeros(nSegments,1);
    segmentPulse = zeros(nSegments,1);
    segmentSamples = zeros(nSegments,segLen);

    rowCounter = 0;

    for iTrial = 1:nTrials
        for iPulse = 1:nPulses

            rowCounter = rowCounter + 1;

            centerIdx = nearestIndex(timeVals,pulseTimes(iPulse));
            idx = (centerIdx + stPoint):(centerIdx + edPoint);

            X(rowCounter,:) = trialData(iTrial,idx);
            segmentTrial(rowCounter) = iTrial;
            segmentPulse(rowCounter) = iPulse;
            segmentSamples(rowCounter,:) = idx;
        end
    end

    D = pairwiseEuclidean(X);

    Xhat = zeros(size(X));

    for iSegment = 1:nSegments

        d = D(iSegment,:);
        [~,order] = sort(d,'ascend');

        % Exclude self
        order(order == iSegment) = [];

        % Exclude all segments from current trial if requested
        if opts.excludeSameTrial
            sameTrial = segmentTrial(order) == segmentTrial(iSegment);
            order(sameTrial) = [];
        end

        if isempty(order)
            order = setdiff(1:nSegments,iSegment);
        end

        Kuse = min(opts.K,length(order));
        neighborIdx = order(1:Kuse);

        templateSegment = median(X(neighborIdx,:),1);

        win = buildSegmentWindow(timeVals,pulseTimes,segmentPulse(iSegment),stPoint,edPoint,segLen);

        Xhat(iSegment,:) = templateSegment .* win;
    end

    artifactEstimate = zeros(nTrials,nSamples);
    cleanedTrials = trialData;

    for iSegment = 1:nSegments
        iTrial = segmentTrial(iSegment);
        idx = segmentSamples(iSegment,:);

        artifactEstimate(iTrial,idx) = artifactEstimate(iTrial,idx) + Xhat(iSegment,:);
        cleanedTrials(iTrial,idx) = cleanedTrials(iTrial,idx) - Xhat(iSegment,:);
    end

    out = struct();
    out.methodName = sprintf('SMARTALiteLocal_K%d_LOTO%d',opts.K,opts.excludeSameTrial);
    out.cleanedData = cleanedTrials;
    out.cleanedTrials = cleanedTrials;
    out.artifactEstimate = artifactEstimate;
    out.trialTemplates = Xhat;
    out.params = opts;
    out.timeVals = timeVals;

    out.diagnostics = struct();
    out.diagnostics.validPulseTimes = pulseTimes;
    out.diagnostics.segmentTrial = segmentTrial;
    out.diagnostics.segmentPulse = segmentPulse;
    out.diagnostics.segmentSamples = segmentSamples;
end

function idx = nearestIndex(t,val)
    [~,idx] = min(abs(t-val));
end

function D = pairwiseEuclidean(X)
    G = sum(X.^2,2);
    D2 = G + G' - 2*(X*X');
    D2(D2 < 0) = 0;
    D = sqrt(D2);
end

function win = buildSegmentWindow(timeVals,pulseTimes,pulseIdx,stPoint,edPoint,segLen)

    win = ones(1,segLen);
    nPulses = length(pulseTimes);

    % Left overlap with previous segment
    if pulseIdx > 1
        prevCenter = nearestIndex(timeVals,pulseTimes(pulseIdx-1));
        thisCenter = nearestIndex(timeVals,pulseTimes(pulseIdx));

        prevEnd = prevCenter + edPoint;
        thisStart = thisCenter + stPoint;

        overlap = prevEnd - thisStart + 1;

        if overlap > 1 && overlap < segLen
            taper = sin(linspace(0,pi/2,overlap)).^2;
            win(1:overlap) = taper;
        end
    end

    % Right overlap with next segment
    if pulseIdx < nPulses
        thisCenter = nearestIndex(timeVals,pulseTimes(pulseIdx));
        nextCenter = nearestIndex(timeVals,pulseTimes(pulseIdx+1));

        thisEnd = thisCenter + edPoint;
        nextStart = nextCenter + stPoint;

        overlap = thisEnd - nextStart + 1;

        if overlap > 1 && overlap < segLen
            taper = cos(linspace(0,pi/2,overlap)).^2;
            win(end-overlap+1:end) = taper;
        end
    end
end