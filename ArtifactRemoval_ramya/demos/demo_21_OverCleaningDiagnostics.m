clear; close all; clc;

%% demo_21_OverCleaningDiagnostics
%
% Goal:
% Diagnose whether cleaned spectra go below the no-stim reference.
%
% Standard FFT error:
%   totalError = norm(cleanedLogFFT - noStimLogFFT)
%
% New diagnostics:
%   aboveReferenceError = norm(max(cleanedLogFFT - noStimLogFFT,0))
%       Residual excess spectral power. Usually leftover artifact.
%
%   belowReferenceError = norm(max(noStimLogFFT - cleanedLogFFT,0))
%       Possible over-cleaning / power loss.
%
%   fractionBelowReference
%       Fraction of frequencies where cleaned spectrum is below no-stim.
%
% Important:
% no-stim is a reference condition, not absolute ground truth.
% But large below-reference error is still a warning sign.

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

% Optional high-frequency diagnostic band
highFreqRange = [80 200];

% PCA settings
artifactWindow = [-0.02 0.4];

% SMARTALite full-cycle settings
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% Over-cleaning weight for balanced score
% lambdaBelow > 1 penalizes going below no-stim more strongly.
lambdaBelow = 1.5;

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
% Required functions
% ========================================================================

requiredMethods = {'PCATemplate','SMARTALite','compute_fft_summary'};

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

fprintf('\nFinal good V1 electrodes used in demo_21:\n');
disp(goodV1Electrodes);
fprintf('Final number of electrodes = %d\n',length(goodV1Electrodes));

%% ========================================================================
% Electrode files
% ========================================================================

elecFiles = dir(fullfile(lfpFolder,'elec*.mat'));

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

methodNames = { ...
    'RawHighStim', ...
    'PCATemplate_K20', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_Ensemble_K3K5'};

methodFieldNames = { ...
    'rawHighStim', ...
    'pcaK20', ...
    'smartaK3', ...
    'smartaK5', ...
    'ensembleK3K5'};

%% ========================================================================
% Output folders
% ========================================================================

resultsFolder = fullfile(folderSourceString,'icms_artifact_removal','results','metrics');

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

figFolder = fullfile(folderSourceString,'icms_artifact_removal','results','figures','demo21_OverCleaningDiagnostics');

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

        %% Run methods
        pcaK20Out = PCATemplate(dataHighStim,timeVals,pcaParamsK20);
        pcaK20Cleaned = getCleanedData(pcaK20Out);

        smartaK3Out = SMARTALite(dataHighStim,timeVals,smartaK3Opts);
        smartaK3Cleaned = getCleanedData(smartaK3Out);

        smartaK5Out = SMARTALite(dataHighStim,timeVals,smartaK5Opts);
        smartaK5Cleaned = getCleanedData(smartaK5Out);

        ensembleCleaned = 0.5 * (smartaK3Cleaned + smartaK5Cleaned);

        dataCell = { ...
            dataHighStim, ...
            pcaK20Cleaned, ...
            smartaK3Cleaned, ...
            smartaK5Cleaned, ...
            ensembleCleaned};

        %% Reference FFT
        fftNoStim = compute_fft_summary(dataNoStim,timeVals,fftIdx);

        freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
                   fftNoStim.freqAxis <= freqRangeForMetric(2);

        highFreqMask = fftNoStim.freqAxis >= highFreqRange(1) & ...
                       fftNoStim.freqAxis <= highFreqRange(2);

        %% Store row
        rowCounter = rowCounter + 1;
        resultRows(rowCounter).electrode = elecNum;

        fprintf('%-28s totalErr  aboveErr  belowErr  fracBelow  balanced\n','Method');

        for iMethod = 1:length(methodNames)

            thisData = dataCell{iMethod};
            thisFFT = compute_fft_summary(thisData,timeVals,fftIdx);

            diagAll = computeOverCleaningDiagnostics(thisFFT,fftNoStim,freqMask,lambdaBelow);
            diagHigh = computeOverCleaningDiagnostics(thisFFT,fftNoStim,highFreqMask,lambdaBelow);

            field = methodFieldNames{iMethod};

            resultRows(rowCounter).([field 'TotalError']) = diagAll.totalError;
            resultRows(rowCounter).([field 'AboveError']) = diagAll.aboveReferenceError;
            resultRows(rowCounter).([field 'BelowError']) = diagAll.belowReferenceError;
            resultRows(rowCounter).([field 'FractionBelow']) = diagAll.fractionBelowReference;
            resultRows(rowCounter).([field 'MeanSignedDiff']) = diagAll.meanSignedDiff;
            resultRows(rowCounter).([field 'BalancedScore']) = diagAll.balancedScore;

            resultRows(rowCounter).([field 'HighTotalError']) = diagHigh.totalError;
            resultRows(rowCounter).([field 'HighAboveError']) = diagHigh.aboveReferenceError;
            resultRows(rowCounter).([field 'HighBelowError']) = diagHigh.belowReferenceError;
            resultRows(rowCounter).([field 'HighFractionBelow']) = diagHigh.fractionBelowReference;
            resultRows(rowCounter).([field 'HighBalancedScore']) = diagHigh.balancedScore;

            fprintf('%-28s %.4f    %.4f    %.4f    %.2f       %.4f\n', ...
                methodNames{iMethod}, ...
                diagAll.totalError, ...
                diagAll.aboveReferenceError, ...
                diagAll.belowReferenceError, ...
                diagAll.fractionBelowReference, ...
                diagAll.balancedScore);
        end

    catch ME

        warning('Failed on elec%d: %s',elecNum,ME.message);

    end
end

if isempty(resultRows)
    error('No electrodes were processed successfully.');
end

resultsTable = struct2table(resultRows);

fprintf('\nDemo 21 over-cleaning diagnostics results table:\n');
disp(resultsTable);

%% ========================================================================
% Summary statistics
% ========================================================================

summaryRows = struct([]);

for iMethod = 1:length(methodNames)

    field = methodFieldNames{iMethod};

    totalError = resultsTable.([field 'TotalError']);
    aboveError = resultsTable.([field 'AboveError']);
    belowError = resultsTable.([field 'BelowError']);
    fractionBelow = resultsTable.([field 'FractionBelow']);
    balancedScore = resultsTable.([field 'BalancedScore']);

    highBelowError = resultsTable.([field 'HighBelowError']);
    highFractionBelow = resultsTable.([field 'HighFractionBelow']);
    highBalancedScore = resultsTable.([field 'HighBalancedScore']);

    summaryRows(iMethod).method = methodNames{iMethod};

    summaryRows(iMethod).meanTotalError = mean(totalError);
    summaryRows(iMethod).stdTotalError = std(totalError);

    summaryRows(iMethod).meanAboveError = mean(aboveError);
    summaryRows(iMethod).stdAboveError = std(aboveError);

    summaryRows(iMethod).meanBelowError = mean(belowError);
    summaryRows(iMethod).stdBelowError = std(belowError);

    summaryRows(iMethod).meanFractionBelow = mean(fractionBelow);

    summaryRows(iMethod).meanBalancedScore = mean(balancedScore);
    summaryRows(iMethod).stdBalancedScore = std(balancedScore);

    summaryRows(iMethod).meanHighBelowError = mean(highBelowError);
    summaryRows(iMethod).meanHighFractionBelow = mean(highFractionBelow);
    summaryRows(iMethod).meanHighBalancedScore = mean(highBalancedScore);
end

summaryTable = struct2table(summaryRows);
summaryTable = sortrows(summaryTable,'meanBalancedScore','ascend');

fprintf('\nDemo 21 over-cleaning diagnostics summary table:\n');
disp(summaryTable);

%% ========================================================================
% Best method counts by different criteria
% ========================================================================

candidateMethodNames = methodNames(2:end);
candidateFieldNames = methodFieldNames(2:end);

totalMatrix = [];
belowMatrix = [];
balancedMatrix = [];

for iMethod = 1:length(candidateFieldNames)
    field = candidateFieldNames{iMethod};

    totalMatrix = [totalMatrix, resultsTable.([field 'TotalError'])];
    belowMatrix = [belowMatrix, resultsTable.([field 'BelowError'])];
    balancedMatrix = [balancedMatrix, resultsTable.([field 'BalancedScore'])];
end

[~,bestTotalIdx] = min(totalMatrix,[],2);
[~,bestBelowIdx] = min(belowMatrix,[],2);
[~,bestBalancedIdx] = min(balancedMatrix,[],2);

bestTotalTable = table( ...
    resultsTable.electrode, ...
    candidateMethodNames(bestTotalIdx)', ...
    candidateMethodNames(bestBelowIdx)', ...
    candidateMethodNames(bestBalancedIdx)', ...
    'VariableNames',{'electrode','bestTotalError','leastBelowReference','bestBalancedScore'} ...
);

fprintf('\nBest method per electrode by criterion:\n');
disp(bestTotalTable);

fprintf('\nBest counts by total FFT error:\n');
printBestCounts(candidateMethodNames,bestTotalIdx);

fprintf('\nBest counts by lowest below-reference error:\n');
printBestCounts(candidateMethodNames,bestBelowIdx);

fprintf('\nBest counts by balanced score:\n');
printBestCounts(candidateMethodNames,bestBalancedIdx);

%% ========================================================================
% Figures
% ========================================================================

% Use candidate methods only, excluding raw
candidateSummary = summaryTable(~strcmp(summaryTable.method,'RawHighStim'),:);

figure('Name','demo21 over-cleaning summary');

barData = [ ...
    candidateSummary.meanAboveError, ...
    candidateSummary.meanBelowError];

bar(categorical(candidateSummary.method),barData);
ylabel('Norm of spectral difference');
title('Above-reference vs below-reference spectral error');
legend({'Above reference: residual excess power','Below reference: possible over-cleaning'}, ...
    'Location','best');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo21_above_vs_below_reference.png'));
    savefig(gcf,fullfile(figFolder,'demo21_above_vs_below_reference.fig'));
end

figure('Name','demo21 balanced score');

bar(categorical(candidateSummary.method),candidateSummary.meanBalancedScore);
ylabel('Balanced score');
title(['Balanced score = aboveError + ' num2str(lambdaBelow) ' * belowError']);
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo21_balanced_score.png'));
    savefig(gcf,fullfile(figFolder,'demo21_balanced_score.fig'));
end

figure('Name','demo21 fraction below reference');

bar(categorical(candidateSummary.method),candidateSummary.meanFractionBelow);
ylabel('Mean fraction of frequencies below no-stim');
title('Fraction of 0-200 Hz spectrum below no-stim reference');
grid on;
xtickangle(30);

if saveFigures
    saveas(gcf,fullfile(figFolder,'demo21_fraction_below_reference.png'));
    savefig(gcf,fullfile(figFolder,'demo21_fraction_below_reference.fig'));
end

%% ========================================================================
% Save results
% ========================================================================

timestampString = datestr(now,'yyyymmdd_HHMMSS');

resultsMatFile = fullfile(resultsFolder,['demo21_overcleaning_diagnostics_results_' timestampString '.mat']);
resultsCsvFile = fullfile(resultsFolder,['demo21_overcleaning_diagnostics_results_' timestampString '.csv']);
summaryCsvFile = fullfile(resultsFolder,['demo21_overcleaning_diagnostics_summary_' timestampString '.csv']);

save(resultsMatFile, ...
    'resultsTable', ...
    'summaryTable', ...
    'bestTotalTable', ...
    'methodNames', ...
    'methodFieldNames', ...
    'goodV1Electrodes', ...
    'pulseTimes', ...
    'lambdaBelow', ...
    'fftWindow', ...
    'freqRangeForMetric', ...
    'highFreqRange');

writetable(resultsTable,resultsCsvFile);
writetable(summaryTable,summaryCsvFile);

fprintf('\nSaved demo_21 results:\n');
fprintf('%s\n',resultsMatFile);
fprintf('%s\n',resultsCsvFile);
fprintf('%s\n',summaryCsvFile);

if saveFigures
    fprintf('\nSaved demo_21 figures in:\n');
    fprintf('%s\n',figFolder);
end

%% ========================================================================
% Final note
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_21 over-cleaning diagnostics complete\n');
fprintf('============================================================\n');
fprintf('Interpretation:\n');
fprintf('Total error tells which method is closest to no-stim overall.\n');
fprintf('Below-reference error warns about possible over-cleaning.\n');
fprintf('Balanced score penalizes below-reference error more strongly.\n');

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

function diagOut = computeOverCleaningDiagnostics(fftCleaned,fftNoStim,freqMask,lambdaBelow)

    diffVals = fftCleaned.logMeanMagnitude(freqMask) - ...
               fftNoStim.logMeanMagnitude(freqMask);

    aboveVals = max(diffVals,0);
    belowVals = max(-diffVals,0);

    diagOut = struct();

    diagOut.totalError = norm(diffVals);
    diagOut.aboveReferenceError = norm(aboveVals);
    diagOut.belowReferenceError = norm(belowVals);

    diagOut.meanAboveReference = mean(aboveVals);
    diagOut.meanBelowReference = mean(belowVals);

    diagOut.fractionBelowReference = mean(diffVals < 0);
    diagOut.meanSignedDiff = mean(diffVals);

    diagOut.balancedScore = diagOut.aboveReferenceError + ...
                            lambdaBelow * diagOut.belowReferenceError;
end

function printBestCounts(methodNames,bestIdx)

    for iMethod = 1:length(methodNames)
        nBest = sum(bestIdx == iMethod);
        fprintf('%s: %d electrodes\n',methodNames{iMethod},nBest);
    end
end