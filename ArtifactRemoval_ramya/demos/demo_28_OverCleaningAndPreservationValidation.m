clear; close all; clc;

%% demo_28_OverCleaningAndPreservationValidation
%
% Purpose:
% Validate artifact-removal methods for both artifact suppression and
% neural-signal preservation.
%
% Validations:
%   1. No-stimulation negative-control distortion
%   2. Pulse-triggered residual artifact around stimulation pulses
%   3. Amplitude-response preservation across stimulation amplitudes
%
% Methods:
%   NoCleaning
%   ERPSubtraction
%   PCATemplate K=10
%   SMARTALite ensemble K3K5
%   SMARTAFull K=3, hp=100 Hz
%
% Main idea:
%   PCA may suppress artifacts strongly but can be aggressive.
%   SMARTAFull/SMARTALite should reduce pulse-locked artifact while
%   preserving broader neural response structure.

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName  = 'dona';
gridType     = 'Microelectrode';
expDate      = '290825';
protocolName = 'GRF_001';

% Base parameter index used in your previous demos:
% noStimCondition   = {1,1,1,5,5,4};
% highStimCondition = {7,1,1,5,5,4};
baseConditionIndex = [1 1 1 5 5 4];

noStimAmpIndex   = 1;
highStimAmpIndex = 7;

% 31 selected V1 electrodes from previous final validation
goodV1Electrodes = [ ...
     3  4  5  6  7  8 10 11 13 ...
    14 15 17 18 19 20 21 22 23 ...
    24 25 26 27 28 29 30 31 34 ...
    35 36 41 42];

% Representative electrode for pulse-triggered and amplitude-tuning plots
representativeElectrode = 7;

% Cleaning/artifact window
artifactWindow = [0 0.4];

% FFT window
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

% ICMS pulse settings
pulseTimes = 0 + (0:7)*0.05;
pulseWindow = [-0.005 0.020];       % pulse-triggered residual window
pulseArtifactMetricWindow = [0 0.006];
pulseRingingMetricWindow  = [0.006 0.020];

% Neural response preservation window
baselineWindow = [-0.20 -0.05];
neuralResponseWindow = [0.45 1.00];

% Methods to test
methodList = { ...
    'NoCleaning', ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_Ensemble_K3K5', ...
    'SMARTAFull_K3_hp100'};

saveFigures = true;

%% ========================================================================
% Paths
% ========================================================================

projectFolder1 = fullfile(folderSourceString,'ICMS_Artifact_Removal');
projectFolder2 = fullfile(folderSourceString,'icms_artifact_removal');

projectMethods = fullfile(folderSourceString,'ICMS_Artifact_Removal','methods');
githubSMARTA = fullfile(folderSourceString,'artifact_removal','Artifact-removal-main');

addpath(githubSMARTA);
addpath(projectMethods,'-begin');

addpath(fullfile(projectFolder1,'methods'));
addpath(fullfile(projectFolder1,'utils'));
% addpath(fullfile(projectFolder1,'external','smarta'));

addpath(fullfile(projectFolder2,'methods'));
addpath(fullfile(projectFolder2,'utils'));
% addpath(fullfile(projectFolder2,'external','smarta'));

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder   = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

resultFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','metrics','demo28_OverCleaningAndPreservation');

figFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','figures','demo28_OverCleaningAndPreservation');

if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

if saveFigures && ~exist(figFolder,'dir')
    mkdir(figFolder);
end

%% ========================================================================
% Load shared files
% ========================================================================

I = load(lfpInfoFile);
P = load(paramFile);

timeVals = I.timeVals(:).';
parameterCombinations = P.parameterCombinations;

if isfield(P,'aValsUnique')
    ampVals = convertParameterValues(P.aValsUnique);
else
    ampVals = 1:size(parameterCombinations,1);
end

noStimIndex = baseConditionIndex;
noStimIndex(1) = noStimAmpIndex;

highStimIndex = baseConditionIndex;
highStimIndex(1) = highStimAmpIndex;

noStimTrials = getTrialsFromParameterCombinations(parameterCombinations,noStimIndex);
highStimTrials = getTrialsFromParameterCombinations(parameterCombinations,highStimIndex);

fprintf('No-stim trials   : %d\n',numel(noStimTrials));
fprintf('High-stim trials : %d\n',numel(highStimTrials));
fprintf('Representative electrode: elec%d\n',representativeElectrode);

fftIdx = find(timeVals > fftWindow(1) & timeVals <= fftWindow(2));

Fs = 1/median(diff(timeVals));
nFFT = length(fftIdx);
freqAxis = (0:nFFT-1) * (Fs/nFFT);

freqMask = freqAxis >= freqRangeForMetric(1) & freqAxis <= freqRangeForMetric(2);

fprintf('FFT samples      : %d\n',nFFT);
fprintf('Frequency step   : %.2f Hz\n',Fs/nFFT);

%% ========================================================================
% Method options
% ========================================================================

opts = struct();
opts.artifactWindow = artifactWindow;
opts.pulseTimes = pulseTimes;
opts.stiFreq = 20;
opts.prePulse = 0.0005;

%% ========================================================================
% 1. No-stimulation negative-control distortion
% ========================================================================

fprintf('\n============================================================\n');
fprintf('1. No-stimulation negative-control distortion\n');
fprintf('============================================================\n');

noStimRows = {};

for iElec = 1:numel(goodV1Electrodes)

    elecNum = goodV1Electrodes(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    if ~exist(elecFile,'file')
        warning('Missing elec%d. Skipping.',elecNum);
        continue;
    end

    fprintf('No-stim distortion: elec%d (%d/%d)\n',elecNum,iElec,numel(goodV1Electrodes));

    D = load(elecFile);

    if ~isfield(D,'analogData')
        warning('analogData missing in elec%d. Skipping.',elecNum);
        continue;
    end

    rawNoStim = D.analogData(noStimTrials,:);
    fftRawNoStim = computeLogFFT(rawNoStim,fftIdx);

    rawTrialSTD = std(rawNoStim,0,1);
    rawMeanTrace = mean(rawNoStim,1);

    for iMethod = 1:numel(methodList)

        methodName = methodList{iMethod};

        cleanedNoStim = applyCleaningMethod(rawNoStim,timeVals,methodName,opts);

        fftCleanNoStim = computeLogFFT(cleanedNoStim,fftIdx);

        diffData = cleanedNoStim - rawNoStim;

        rmsDistortion = safeRMS(diffData(:));
        normalizedRMSDistortion = rmsDistortion / max(safeRMS(rawNoStim(:)),eps);

        fftDistortion = norm(fftCleanNoStim(freqMask) - fftRawNoStim(freqMask));

        cleanTrialSTD = std(cleanedNoStim,0,1);
        trialSTDDistortion = safeRMS(cleanTrialSTD - rawTrialSTD);

        cleanMeanTrace = mean(cleanedNoStim,1);
        meanTraceCorrelation = safeCorr(rawMeanTrace(:),cleanMeanTrace(:));

        noStimRows(end+1,:) = { ...
            elecNum, ...
            methodName, ...
            rmsDistortion, ...
            normalizedRMSDistortion, ...
            fftDistortion, ...
            trialSTDDistortion, ...
            meanTraceCorrelation};
    end
end

noStimTable = cell2table(noStimRows, ...
    'VariableNames',{ ...
    'electrode', ...
    'method', ...
    'rmsDistortion', ...
    'normalizedRMSDistortion', ...
    'fftDistortion', ...
    'trialSTDDistortion', ...
    'meanTraceCorrelation'});

noStimSummary = summarizeByMethod(noStimTable,methodList);

disp('No-stim distortion summary:');
disp(noStimSummary);

writetable(noStimTable,fullfile(resultFolder,'demo28_no_stim_distortion_long.csv'));
writetable(noStimSummary,fullfile(resultFolder,'demo28_no_stim_distortion_summary.csv'));

%% Plot no-stim distortion summary

fig1 = figure('Color','w','Position',[100 100 1000 500]);

bar(noStimSummary.meanNormalizedRMSDistortion);
set(gca,'XTick',1:numel(methodList),'XTickLabel',prettyMethodNames(methodList));
xtickangle(25);
ylabel('Mean normalized RMS distortion');
title('No-stimulation negative-control distortion');
grid on;
box off;

if saveFigures
    saveAllFigureFormats(fig1,figFolder,'demo28_no_stim_negative_control_distortion');
end

%% ========================================================================
% 2. Pulse-triggered residual artifact on representative electrode
% ========================================================================

fprintf('\n============================================================\n');
fprintf('2. Pulse-triggered residual artifact: elec%d\n',representativeElectrode);
fprintf('============================================================\n');

repFile = fullfile(lfpFolder,['elec' num2str(representativeElectrode) '.mat']);
Drep = load(repFile);

dataNoStimRep = Drep.analogData(noStimTrials,:);
dataHighRep   = Drep.analogData(highStimTrials,:);

pulseRows = {};
pulseTraces = struct();

for iMethod = 1:numel(methodList)

    methodName = methodList{iMethod};

    fprintf('Pulse residual method: %s\n',methodName);

    if strcmp(methodName,'NoCleaning')
        cleanedHigh = dataHighRep;
        displayName = 'Raw high-stim';
    else
        cleanedHigh = applyCleaningMethod(dataHighRep,timeVals,methodName,opts);
        displayName = prettyMethodName(methodName);
    end

    [relT,avgPulseTrace] = computePulseTriggeredAverage( ...
        cleanedHigh,timeVals,pulseTimes,pulseWindow);

    artIdx = relT >= pulseArtifactMetricWindow(1) & relT <= pulseArtifactMetricWindow(2);
    ringIdx = relT >= pulseRingingMetricWindow(1) & relT <= pulseRingingMetricWindow(2);

    pulseArtifactRMS = safeRMS(avgPulseTrace(artIdx));
    pulseArtifactPeakToPeak = max(avgPulseTrace(artIdx)) - min(avgPulseTrace(artIdx));

    pulseRingingRMS = safeRMS(avgPulseTrace(ringIdx));
    pulseOverallRMS = safeRMS(avgPulseTrace);

    pulseRows(end+1,:) = { ...
        methodName, ...
        pulseArtifactRMS, ...
        pulseArtifactPeakToPeak, ...
        pulseRingingRMS, ...
        pulseOverallRMS};

    pulseTraces(iMethod).method = methodName;
    pulseTraces(iMethod).displayName = displayName;
    pulseTraces(iMethod).relT = relT;
    pulseTraces(iMethod).avgPulseTrace = avgPulseTrace;
end

pulseTable = cell2table(pulseRows, ...
    'VariableNames',{ ...
    'method', ...
    'pulseArtifactRMS_0to6ms', ...
    'pulseArtifactPeakToPeak_0to6ms', ...
    'pulseRingingRMS_6to20ms', ...
    'pulseOverallRMS'});

disp('Pulse-triggered residual summary:');
disp(pulseTable);

writetable(pulseTable,fullfile(resultFolder, ...
    ['demo28_pulse_triggered_residual_elec' num2str(representativeElectrode) '.csv']));

%% Plot pulse-triggered residual

fig2 = figure('Color','w','Position',[100 100 1000 550]);

hold on;

for iMethod = 1:numel(pulseTraces)
    plot(1000*pulseTraces(iMethod).relT, ...
        pulseTraces(iMethod).avgPulseTrace, ...
        'LineWidth',1.6);
end

xline(0,'--','Pulse','LineWidth',1.0);

xlabel('Time from pulse (ms)');
ylabel('Pulse-triggered mean LFP (a.u.)');
title(['Pulse-triggered residual artifact: elec' num2str(representativeElectrode)]);
legend({pulseTraces.displayName},'Location','best');
grid on;
box off;

if saveFigures
    saveAllFigureFormats(fig2,figFolder, ...
        ['demo28_pulse_triggered_residual_elec' num2str(representativeElectrode)]);
end

%% ========================================================================
% 3. Amplitude-response preservation on representative electrode
% ========================================================================

fprintf('\n============================================================\n');
fprintf('3. Amplitude-response preservation: elec%d\n',representativeElectrode);
fprintf('============================================================\n');

ampRows = {};

nAmp = size(parameterCombinations,1);

for iAmp = 1:nAmp

    thisIndex = baseConditionIndex;
    thisIndex(1) = iAmp;

    thisTrials = getTrialsFromParameterCombinations(parameterCombinations,thisIndex);

    if isempty(thisTrials)
        continue;
    end

    if iAmp <= numel(ampVals)
        thisAmpVal = ampVals(iAmp);
    else
        thisAmpVal = iAmp;
    end

    dataThisAmp = Drep.analogData(thisTrials,:);

    for iMethod = 1:numel(methodList)

        methodName = methodList{iMethod};

        % Do not clean 0 uA/no-stim condition for preservation curve.
        % There is no stimulation artifact there.
        if thisAmpVal == 0 || strcmp(methodName,'NoCleaning')
            cleanedData = dataThisAmp;
        else
            cleanedData = applyCleaningMethod(dataThisAmp,timeVals,methodName,opts);
        end

        [responseRMS,responseMean,responseAUC] = computeNeuralResponseMetric( ...
            cleanedData,timeVals,baselineWindow,neuralResponseWindow);

        ampRows(end+1,:) = { ...
            iAmp, ...
            thisAmpVal, ...
            methodName, ...
            responseRMS, ...
            responseMean, ...
            responseAUC};
    end
end

ampTable = cell2table(ampRows, ...
    'VariableNames',{ ...
    'ampIndex', ...
    'amplitudeMicroAmps', ...
    'method', ...
    'responseRMS', ...
    'responseMean', ...
    'responseAUC'});

disp('Amplitude-response preservation table:');
disp(ampTable);

writetable(ampTable,fullfile(resultFolder, ...
    ['demo28_amplitude_response_preservation_elec' num2str(representativeElectrode) '.csv']));

%% Plot amplitude response

fig3 = figure('Color','w','Position',[100 100 1000 550]);

hold on;

for iMethod = 1:numel(methodList)

    methodName = methodList{iMethod};
    idx = strcmp(ampTable.method,methodName);

    x = ampTable.amplitudeMicroAmps(idx);
    y = ampTable.responseRMS(idx);

    [xSort,order] = sort(x);
    ySort = y(order);

    plot(xSort,ySort,'-o','LineWidth',1.6,'MarkerSize',6);
end

xlabel('Stimulation amplitude (\muA)');
ylabel('Post-stimulation response RMS (a.u.)');
title(['Amplitude-response preservation: elec' num2str(representativeElectrode)]);
legend(prettyMethodNames(methodList),'Location','best');
grid on;
box off;

if saveFigures
    saveAllFigureFormats(fig3,figFolder, ...
        ['demo28_amplitude_response_preservation_elec' num2str(representativeElectrode)]);
end

%% ========================================================================
% Final summary
% ========================================================================

fprintf('\n============================================================\n');
fprintf('demo_28 complete\n');
fprintf('============================================================\n');
fprintf('Saved tables in:\n%s\n',resultFolder);
fprintf('Saved figures in:\n%s\n',figFolder);

fprintf('\nRecommended interpretation:\n');
fprintf('1. Low no-stim distortion supports signal preservation.\n');
fprintf('2. Low pulse-triggered residual supports artifact suppression.\n');
fprintf('3. Preserved amplitude-response tuning supports physiological preservation.\n');

%% ========================================================================
% Helper functions
% ========================================================================

function trials = getTrialsFromParameterCombinations(parameterCombinations,indexVector)

    dims = size(parameterCombinations);
    nDims = ndims(parameterCombinations);

    if numel(indexVector) < nDims
        indexVector = [indexVector ones(1,nDims-numel(indexVector))];
    end

    indexVector = indexVector(1:nDims);

    for i = 1:nDims
        if indexVector(i) > dims(i)
            error('Index %d exceeds parameterCombinations dimension %d.',indexVector(i),i);
        end
    end

    indexCell = num2cell(indexVector);
    trials = parameterCombinations{indexCell{:}};
end

function vals = convertParameterValues(valsIn)

    vals = valsIn;

    vals(vals > 16384) = vals(vals > 16384) - 32768;
end

function signalClean = applyCleaningMethod(signalRaw,timeVals,methodName,opts)

    switch methodName

        case 'NoCleaning'
            signalClean = signalRaw;

        case 'ERPSubtraction'
            erpParams = struct();
            erpParams.subtractWindow = opts.artifactWindow;
            erpParams.doBaselineCorrection = false;

            out = ERPSubtraction(signalRaw,timeVals,erpParams);
            signalClean = getCleanedData(out);

        case 'PCATemplate_K10'
            pcaParams = struct();
            pcaParams.artifactWindow = opts.artifactWindow;
            pcaParams.numComponents = 10;
            pcaParams.removeMeanTemplate = true;
            pcaParams.taperEdgeMS = 2;
            pcaParams.doBaselineCorrection = false;

            out = PCATemplate(signalRaw,timeVals,pcaParams);
            signalClean = getCleanedData(out);

        case 'SMARTALite_Ensemble_K3K5'
            liteOpts = struct();
            liteOpts.pulseTimes = opts.pulseTimes;
            liteOpts.stiFreq = opts.stiFreq;
            liteOpts.prePulse = opts.prePulse;
            liteOpts.window = opts.artifactWindow;
            liteOpts.computeFFT = false;
            liteOpts.excludeSameTrial = true;

            liteOpts.K = 3;
            outK3 = SMARTALite(signalRaw,timeVals,liteOpts);
            cleanK3 = getCleanedData(outK3);

            liteOpts.K = 5;
            outK5 = SMARTALite(signalRaw,timeVals,liteOpts);
            cleanK5 = getCleanedData(outK5);

            signalClean = 0.5 * (cleanK3 + cleanK5);

        case 'SMARTAFull_K3_hp100'
            fullOpts = struct();
            fullOpts.pulseTimes = opts.pulseTimes;
            fullOpts.stiFreq = opts.stiFreq;
            fullOpts.prePulse = opts.prePulse;
            fullOpts.window = opts.artifactWindow;
            fullOpts.computeFFT = false;
            fullOpts.excludeSameTrial = true;

            fullOpts.K = 3;
            fullOpts.hpCutoff = 100;
            fullOpts.shrinkLoss = 'fro';
            fullOpts.shrinkKL = 10;
            fullOpts.shrinkKH = 15;

            out = SMARTAFull(signalRaw,timeVals,fullOpts);
            signalClean = getCleanedData(out);

        otherwise
            error('Unknown methodName: %s',methodName);
    end
end

function cleanedData = getCleanedData(methodOut)

    if isfield(methodOut,'cleanedData')
        cleanedData = methodOut.cleanedData;
    elseif isfield(methodOut,'cleanedTrials')
        cleanedData = methodOut.cleanedTrials;
    elseif isfield(methodOut,'z')
        cleanedData = methodOut.z;
    else
        error('Method output has no cleaned data field.');
    end
end

function y = computeLogFFT(data,fftIdx)

    fftVals = abs(fft(data(:,fftIdx),[],2));
    y = log10(mean(fftVals,1) + eps);
end

function r = safeRMS(x)

    x = x(:);
    r = sqrt(mean(x.^2,'omitnan'));
end

function c = safeCorr(x,y)

    x = x(:);
    y = y(:);

    good = isfinite(x) & isfinite(y);

    if sum(good) < 3
        c = NaN;
        return;
    end

    x = x(good);
    y = y(good);

    if std(x) == 0 || std(y) == 0
        c = NaN;
        return;
    end

    C = corrcoef(x,y);
    c = C(1,2);
end

function summaryTable = summarizeByMethod(T,methodList)

    rows = {};

    for iMethod = 1:numel(methodList)

        methodName = methodList{iMethod};
        idx = strcmp(T.method,methodName);

        rows(end+1,:) = { ...
            methodName, ...
            sum(idx), ...
            mean(T.rmsDistortion(idx),'omitnan'), ...
            std(T.rmsDistortion(idx),'omitnan'), ...
            mean(T.normalizedRMSDistortion(idx),'omitnan'), ...
            std(T.normalizedRMSDistortion(idx),'omitnan'), ...
            mean(T.fftDistortion(idx),'omitnan'), ...
            std(T.fftDistortion(idx),'omitnan'), ...
            mean(T.trialSTDDistortion(idx),'omitnan'), ...
            std(T.trialSTDDistortion(idx),'omitnan'), ...
            mean(T.meanTraceCorrelation(idx),'omitnan')};
    end

    summaryTable = cell2table(rows, ...
        'VariableNames',{ ...
        'method', ...
        'nElectrodes', ...
        'meanRMSDistortion', ...
        'stdRMSDistortion', ...
        'meanNormalizedRMSDistortion', ...
        'stdNormalizedRMSDistortion', ...
        'meanFFTDistortion', ...
        'stdFFTDistortion', ...
        'meanTrialSTDDistortion', ...
        'stdTrialSTDDistortion', ...
        'meanTraceCorrelation'});
end

function [relT,avgPulseTrace] = computePulseTriggeredAverage(data,timeVals,pulseTimes,pulseWindow)

    Fs = 1/median(diff(timeVals));

    relT = pulseWindow(1):1/Fs:pulseWindow(2);
    nSamp = numel(relT);

    allSegments = [];

    for iPulse = 1:numel(pulseTimes)

        centerTime = pulseTimes(iPulse);
        sampleTimes = centerTime + relT;

        if sampleTimes(1) < timeVals(1) || sampleTimes(end) > timeVals(end)
            continue;
        end

        seg = interp1(timeVals,data.',sampleTimes,'linear','extrap').';
        allSegments = cat(1,allSegments,seg);
    end

    if isempty(allSegments)
        avgPulseTrace = nan(1,nSamp);
        return;
    end

    baselineIdx = relT < -0.001;

    if any(baselineIdx)
        allSegments = allSegments - mean(allSegments(:,baselineIdx),2);
    end

    avgPulseTrace = mean(allSegments,1,'omitnan');
end

function [responseRMS,responseMean,responseAUC] = computeNeuralResponseMetric(data,timeVals,baselineWindow,responseWindow)

    baselineIdx = timeVals >= baselineWindow(1) & timeVals <= baselineWindow(2);
    responseIdx = timeVals >= responseWindow(1) & timeVals <= responseWindow(2);

    meanTrace = mean(data,1,'omitnan');

    baselineVal = mean(meanTrace(baselineIdx),'omitnan');
    responseTrace = meanTrace(responseIdx) - baselineVal;

    responseRMS = safeRMS(responseTrace);
    responseMean = mean(responseTrace,'omitnan');

    dt = median(diff(timeVals));
    responseAUC = sum(abs(responseTrace),'omitnan') * dt;
end

function names = prettyMethodNames(methodList)

    names = cell(size(methodList));

    for i = 1:numel(methodList)
        names{i} = prettyMethodName(methodList{i});
    end
end

function name = prettyMethodName(methodName)

    switch methodName
        case 'NoCleaning'
            name = 'No cleaning / raw';
        case 'ERPSubtraction'
            name = 'ERPSubtraction';
        case 'PCATemplate_K10'
            name = 'PCATemplate K=10';
        case 'SMARTALite_Ensemble_K3K5'
            name = 'SMARTALite ensemble K3K5';
        case 'SMARTAFull_K3_hp100'
            name = 'SMARTAFull K=3, hp=100 Hz';
        otherwise
            name = methodName;
    end
end

function saveAllFigureFormats(figHandle,figFolder,baseName)

    pngFile = fullfile(figFolder,[baseName '.png']);
    figFile = fullfile(figFolder,[baseName '.fig']);
    pdfFile = fullfile(figFolder,[baseName '.pdf']);

    saveas(figHandle,pngFile);
    savefig(figHandle,figFile);

    try
        exportgraphics(figHandle,pdfFile,'ContentType','vector');
    catch
        exportgraphics(figHandle,pdfFile,'ContentType','image');
    end

    fprintf('Saved figure:\n%s\n%s\n%s\n',pngFile,figFile,pdfFile);
end