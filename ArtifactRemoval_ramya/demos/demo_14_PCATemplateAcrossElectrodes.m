clear; close all; clc;

%% demo_14_PCATemplateAcrossElectrodes
%
% Goal:
% Test whether PCATemplate works beyond elec1.
%
% Corrected electrode selection:
% 1. Load highRMSElectrodes from donaMicroelectrodeRFData.mat
% 2. Optionally remove bad impedance electrodes, if impedanceValues.mat exists
% 3. Restrict to V1 electrodes 1:48
% 4. Exclude stimulation electrode elec1
%
% Expected for your current data:
% highRMSElectrodes across full array = 66 electrodes
% good V1 electrodes before stim exclusion = 32 electrodes
% final electrodes after excluding elec1 = 31 electrodes

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

% PCA settings
artifactWindow = [-0.02 0.4];
pcaKConservative = 3;
pcaKAggressive = 20;

% FFT settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% Run settings
% Use Inf for all electrodes.
% Use smaller number, e.g. 10, for quick testing.
maxElectrodesToRun = Inf;

% Electrode selection
v1Electrodes = 1:48;
stimElectrode = 1;
excludeStimElectrode = true;

%% ========================================================================
% Paths
% ========================================================================

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder   = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

I = load(lfpInfoFile);
P = load(paramFile);

timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

%% ========================================================================
% Get good electrodes
% ========================================================================

% The good electrode selection:
% badChannels = impedance channels above cutoff or NaN
% goodElectrodes = setdiff(highRMSElectrodes,badChannels)
%
% We additionally restrict to V1 electrodes 1:48.

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

% Same good-electrode definition as the project analysis pipeline
goodElectrodesAll = setdiff(highRMSElectrodes,badChannels);

% Restrict to V1
goodV1Electrodes = intersect(goodElectrodesAll,v1Electrodes);

fprintf('\nRFData file used:\n');
disp(rfDataPath);

fprintf('\nHigh-RMS electrodes from RFData:\n');
disp(highRMSElectrodes);

fprintf('\nBad impedance channels:\n');
disp(badChannels);

fprintf('\nGood electrodes across full array:\n');
disp(goodElectrodesAll);
fprintf('Number of good electrodes across full array = %d\n',length(goodElectrodesAll));

fprintf('\nGood V1 electrodes before stim-electrode exclusion:\n');
disp(goodV1Electrodes);
fprintf('Number of good V1 electrodes = %d\n',length(goodV1Electrodes));

if excludeStimElectrode
    goodV1Electrodes = setdiff(goodV1Electrodes,stimElectrode);
    fprintf('\nExcluding stimulation electrode elec%d.\n',stimElectrode);
end

fprintf('\nFinal good V1 electrodes used in demo_14:\n');
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

    % Anchored regex avoids hidden files like ._elec1.mat
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

%% Restrict files to corrected good V1 electrode list

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
% Trial indices and FFT indices
% ========================================================================

noStimTrials = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

fftIdx = find(timeVals > fftWindow(1) & timeVals < fftWindow(2));

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));
fprintf('FFT window has %d samples.\n',length(fftIdx));
fprintf('Window duration = %.4f s.\n',timeVals(fftIdx(end))-timeVals(fftIdx(1)));

%% ========================================================================
% Method parameters
% ========================================================================

% ERPSubtraction
erpParams = struct();
erpParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
erpParams.doBaselineCorrection = false;

% ERPAligned
alignedParams = struct();
alignedParams.subtractWindow = [timeVals(1) timeVals(end)+eps];
alignedParams.alignWindow = [-0.01 0.03];
alignedParams.maxShiftMS = 10;
alignedParams.doBaselineCorrection = false;

% Tuned hybrid
hybridPulseParams = struct();
hybridPulseParams.pulseSearchWindow = [-0.02 0.32];
hybridPulseParams.pulseWindowMS = [-10 40];
hybridPulseParams.minPulseDistanceMS = 25;
hybridPulseParams.expectedNumPulses = 7;
hybridPulseParams.maxNumPulses = 20;
hybridPulseParams.thresholdMAD = 2;
hybridPulseParams.localBaselineMS = 20;
hybridPulseParams.templateStatistic = 'median';
hybridPulseParams.taperEdgeMS = 2;
hybridPulseParams.doBaselineCorrection = false;

hybridParams = struct();
hybridParams.alignedParams = alignedParams;
hybridParams.pulseParams = hybridPulseParams;

% PCA K=3
pcaParamsK3 = struct();
pcaParamsK3.artifactWindow = artifactWindow;
pcaParamsK3.numComponents = pcaKConservative;
pcaParamsK3.removeMeanTemplate = true;
pcaParamsK3.taperEdgeMS = 2;
pcaParamsK3.doBaselineCorrection = false;

% PCA K=20
pcaParamsK20 = struct();
pcaParamsK20.artifactWindow = artifactWindow;
pcaParamsK20.numComponents = pcaKAggressive;
pcaParamsK20.removeMeanTemplate = true;
pcaParamsK20.taperEdgeMS = 2;
pcaParamsK20.doBaselineCorrection = false;

%% ========================================================================
% Loop across electrodes
% ========================================================================

resultRows = struct([]);
rowCounter = 0;

for iElec = 1:length(elecFiles)

    elecNum = elecNums(iElec);
    elecFile = fullfile(lfpFolder,elecFiles(iElec).name);

    fprintf('\nProcessing electrode %d/%d: elec%d\n',iElec,length(elecFiles),elecNum);

    try
        D = load(elecFile);

        if ~isfield(D,'analogData')
            warning('Skipping elec%d: analogData not found.',elecNum);
            continue;
        end

        analogData = D.analogData;

        dataNoStim = analogData(noStimTrials,:);
        dataHighStim = analogData(highStimTrials,:);

        %% FFT reference for this electrode
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

        %% ERPAligned
        alignedOut = ERPAligned(dataHighStim,timeVals,alignedParams);
        fftAligned = compute_fft_summary(alignedOut.cleanedData,timeVals,fftIdx);

        alignedError = norm(fftAligned.logMeanMagnitude(freqMask) - ...
                            fftNoStim.logMeanMagnitude(freqMask));

        alignedImprovement = 100 * (rawError - alignedError) / rawError;

        %% Tuned hybrid
        hybridOut = ERPAlignedPulsewise(dataHighStim,timeVals,hybridParams);
        fftHybrid = compute_fft_summary(hybridOut.cleanedData,timeVals,fftIdx);

        hybridError = norm(fftHybrid.logMeanMagnitude(freqMask) - ...
                           fftNoStim.logMeanMagnitude(freqMask));

        hybridImprovement = 100 * (rawError - hybridError) / rawError;

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

        %% Store
        rowCounter = rowCounter + 1;

        resultRows(rowCounter).electrode = elecNum;
        resultRows(rowCounter).rawError = rawError;

        resultRows(rowCounter).erpError = erpError;
        resultRows(rowCounter).erpImprovement = erpImprovement;

        resultRows(rowCounter).alignedError = alignedError;
        resultRows(rowCounter).alignedImprovement = alignedImprovement;

        resultRows(rowCounter).hybridError = hybridError;
        resultRows(rowCounter).hybridImprovement = hybridImprovement;

        resultRows(rowCounter).pcaK3Error = pcaK3Error;
        resultRows(rowCounter).pcaK3Improvement = pcaK3Improvement;

        resultRows(rowCounter).pcaK20Error = pcaK20Error;
        resultRows(rowCounter).pcaK20Improvement = pcaK20Improvement;

        resultRows(rowCounter).alignedClippedTrials = sum(alignedOut.wasClipped);
        resultRows(rowCounter).hybridDetectedPulses = length(hybridOut.pulseOut.pulseTimes);

        fprintf('Raw error       = %.4f\n',rawError);
        fprintf('ERP             = %.4f, improvement %.2f %%\n',erpError,erpImprovement);
        fprintf('ERPAligned      = %.4f, improvement %.2f %%\n',alignedError,alignedImprovement);
        fprintf('Hybrid          = %.4f, improvement %.2f %%\n',hybridError,hybridImprovement);
        fprintf('PCA K=3         = %.4f, improvement %.2f %%\n',pcaK3Error,pcaK3Improvement);
        fprintf('PCA K=20        = %.4f, improvement %.2f %%\n',pcaK20Error,pcaK20Improvement);

    catch ME

        warning('Failed on elec%d: %s',elecNum,ME.message);

    end
end

if isempty(resultRows)
    error('No electrodes were processed successfully.');
end

resultsTable = struct2table(resultRows);

fprintf('\nAcross-electrode results table:\n');
disp(resultsTable);

%% ========================================================================
% Summary statistics
% ========================================================================

methodNames = { ...
    'ERPSubtraction', ...
    'ERPAligned', ...
    'TunedHybrid', ...
    'PCATemplate_K3', ...
    'PCATemplate_K20'};

improvementMatrix = [ ...
    resultsTable.erpImprovement, ...
    resultsTable.alignedImprovement, ...
    resultsTable.hybridImprovement, ...
    resultsTable.pcaK3Improvement, ...
    resultsTable.pcaK20Improvement];

errorMatrix = [ ...
    resultsTable.erpError, ...
    resultsTable.alignedError, ...
    resultsTable.hybridError, ...
    resultsTable.pcaK3Error, ...
    resultsTable.pcaK20Error];

meanImprovement = mean(improvementMatrix,1);
stdImprovement = std(improvementMatrix,0,1);

meanError = mean(errorMatrix,1);
stdError = std(errorMatrix,0,1);

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

fprintf('\nAcross-electrode summary table:\n');
disp(summaryTable);

%% ========================================================================
% Count best method per electrode
% ========================================================================

[~,bestMethodIdx] = min(errorMatrix,[],2);

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
% Figure 1: Mean improvement across electrodes
% ========================================================================

figure('Name','Across electrodes: Mean improvement');

bar(meanImprovement);
hold on;
errorbar(1:length(methodNames),meanImprovement,stdImprovement,'k.','LineWidth',1.2);

xticks(1:length(methodNames));
xticklabels(methodNames);
xtickangle(30);

ylabel('Improvement relative to raw high-stim (%)');
title('Mean artifact-removal improvement across good V1 electrodes');
grid on;

%% ========================================================================
% Figure 2: Error distributions across electrodes
% ========================================================================

figure('Name','Across electrodes: Error distributions');

boxplot(errorMatrix,'Labels',methodNames);
ylabel('FFT error vs no-stim');
title('Error distribution across good V1 electrodes');
xtickangle(30);
grid on;

%% ========================================================================
% Figure 3: Improvement distributions across electrodes
% ========================================================================

figure('Name','Across electrodes: Improvement distributions');

boxplot(improvementMatrix,'Labels',methodNames);
ylabel('Improvement relative to raw high-stim (%)');
title('Improvement distribution across good V1 electrodes');
xtickangle(30);
grid on;

%% ========================================================================
% Figure 4: Electrode-wise improvements
% ========================================================================

figure('Name','Across electrodes: Electrode-wise improvements');

plot(resultsTable.electrode,resultsTable.erpImprovement,'o-','LineWidth',1.1);
hold on;
plot(resultsTable.electrode,resultsTable.alignedImprovement,'o-','LineWidth',1.1);
plot(resultsTable.electrode,resultsTable.hybridImprovement,'o-','LineWidth',1.1);
plot(resultsTable.electrode,resultsTable.pcaK3Improvement,'o-','LineWidth',1.1);
plot(resultsTable.electrode,resultsTable.pcaK20Improvement,'o-','LineWidth',1.1);

xlabel('Electrode');
ylabel('Improvement (%)');
title('Artifact-removal improvement by good V1 electrode');
legend(methodNames,'Location','best');
grid on;

%% ========================================================================
% Figure 5: PCA K=20 vs PCA K=3
% ========================================================================

figure('Name','Across electrodes: PCA K20 vs K3');

scatter(resultsTable.pcaK3Improvement,resultsTable.pcaK20Improvement,60,resultsTable.electrode,'filled');
colorbar;
xlabel('PCATemplate K=3 improvement (%)');
ylabel('PCATemplate K=20 improvement (%)');
title('PCA K=20 vs K=3 across good V1 electrodes');
grid on;

hold on;
minVal = min([resultsTable.pcaK3Improvement; resultsTable.pcaK20Improvement]);
maxVal = max([resultsTable.pcaK3Improvement; resultsTable.pcaK20Improvement]);
plot([minVal maxVal],[minVal maxVal],'k--','LineWidth',1.2);

%% ========================================================================
% Save results
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo14_goodV1_across_electrodes_results_' timestampString '.mat']);
resultsCsvFile = fullfile(resultsFolder,['demo14_goodV1_across_electrodes_results_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo14_goodV1_across_electrodes_summary_' timestampString '.csv']);

save(resultsMatFile, ...
    'resultsTable', ...
    'summaryTable', ...
    'bestMethodTable', ...
    'methodNames', ...
    'goodV1Electrodes', ...
    'goodElectrodesAll', ...
    'highRMSElectrodes', ...
    'badChannels', ...
    'rfDataPath');

writetable(resultsTable,resultsCsvFile);
writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',resultsCsvFile);
fprintf('%s\n',summaryCsvFile);

%% ========================================================================
% Final printed interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('Across-electrode PCATemplate summary: corrected good V1 list\n');
fprintf('============================================================\n');

fprintf('Number of electrodes processed = %d\n',height(resultsTable));
fprintf('Artifact window = [%.3f %.3f] s\n',artifactWindow(1),artifactWindow(2));
fprintf('FFT metric range = [%d %d] Hz\n',freqRangeForMetric(1),freqRangeForMetric(2));

fprintf('\nElectrodes processed:\n');
disp(resultsTable.electrode');

fprintf('\nMethod ranking by mean FFT error:\n');
disp(summaryTable);

fprintf('\nInterpretation:\n');
fprintf('This corrected demo_14 uses high-RMS V1 electrodes from RFData, excluding elec1.\n');
fprintf('If PCATemplate K=20 has the lowest mean error and wins on most electrodes,\n');
fprintf('then the PCA result is stronger than the single-electrode result.\n');
fprintf('If PCATemplate K=3 is close to K=20, K=3 may be the safer conservative choice.\n');