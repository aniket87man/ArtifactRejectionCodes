clear; close all; clc;

%% demo_27_BestSMARTAFFTRepresentative
%
% Purpose:
% Find and plot a strong representative FFT example where SMARTALite or
% SMARTAFull gives the best-looking spectral cleaning result.
%
% Important:
% This is a representative/development figure, not the final 31-electrode
% validation. By default it searches elec1 because the earlier strong SMARTA
% plot was generated from elec1, the stimulation electrode.
%
% Main output:
%   1. Best SMARTA/SMARTAFull FFT comparison figure
%   2. SMARTALite vs SMARTAFull comparison figure
%   3. Search table of FFT/harmonic errors

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

projectMethods = fullfile(folderSourceString,'ICMS_Artifact_Removal','methods');
githubSMARTA = fullfile(folderSourceString,'artifact_removal','Artifact-removal-main');

addpath(githubSMARTA);
addpath(projectMethods,'-begin');

subjectName  = 'dona';
gridType     = 'Microelectrode';
expDate      = '290825';
protocolName = 'GRF_001';

% Conditions
noStimCondition   = {1,1,1,5,5,4};
highStimCondition = {7,1,1,5,5,4};

% ------------------------------------------------------------------------
% For reproducing the earlier strong-looking SMARTA result, use elec1.
% If you want to search more electrodes later, replace with:
% electrodesToSearch = [1 3 4 5 6 7 8 10 11 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 34 35 36 41 42];
% ------------------------------------------------------------------------
electrodesToSearch = 7;

% FFT settings
fftWindow = [0 0.4];
freqRangeForPlot = [0 200];
freqRangeForMetric = [0 200];

% Harmonic metric settings
stimFreq = 20;
harmonicsToUse = stimFreq:stimFreq:200;
harmonicHalfBandwidthHz = 1.0;

% Artifact-removal window
artifactWindow = [0 0.4];

% Pulse settings
pulseTimes = 0 + (0:7)*0.05;
smartaStimFreq = 20;
smartaPrePulse = 0.0005;

% SMARTALite search
liteKList = [3 5 10 20 40 60];

% SMARTAFull search
fullKList = [3 5 10 15 20];
fullHpCutoffs = [100 150 200 300];

% Search objective
% Larger harmonicWeight favors plots with fewer visible stimulation harmonics.
harmonicWeight = 3.0;

% Whether to include PCA in the final plot.
% Keep false if the goal is a SMARTA-focused representative figure.
showPCAInFinalFigure = false;

saveFigures = true;

%% ========================================================================
% Paths
% ========================================================================

projectFolder = fullfile(folderSourceString,'ICMS_Artifact_Removal');
altProjectFolder = fullfile(folderSourceString,'icms_artifact_removal');

addpath(fullfile(projectFolder,'methods'));
addpath(fullfile(projectFolder,'utils'));
addpath(fullfile(projectFolder,'external','smarta'));

addpath(fullfile(altProjectFolder,'methods'));
addpath(fullfile(altProjectFolder,'utils'));
addpath(fullfile(altProjectFolder,'external','smarta'));

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFolder   = fullfile(baseFolder,'segmentedData','lfp');
lfpInfoFile = fullfile(lfpFolder,'lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

figFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','figures','demo27_BestSMARTAFullFFTRepresentative');

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

noStimTrials = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

fftIdx = find(timeVals > fftWindow(1) & timeVals <= fftWindow(2));
Fs = 1/median(diff(timeVals));
nFFT = length(fftIdx);
freqAxis = (0:nFFT-1) * (Fs/nFFT);

freqMask = freqAxis >= freqRangeForMetric(1) & freqAxis <= freqRangeForMetric(2);
plotMask = freqAxis >= freqRangeForPlot(1) & freqAxis <= freqRangeForPlot(2);

harmonicMask = buildHarmonicMask(freqAxis,harmonicsToUse,harmonicHalfBandwidthHz,freqRangeForMetric);

fprintf('No-stim trials   : %d\n',length(noStimTrials));
fprintf('High-stim trials : %d\n',length(highStimTrials));
fprintf('FFT samples      : %d\n',length(fftIdx));
fprintf('Frequency step   : %.2f Hz\n',Fs/nFFT);

%% ========================================================================
% Check functions
% ========================================================================

requiredMethods = {'ERPSubtraction','SMARTALite','SMARTAFull'};

optionalMethods = {'TrialWiseTemplateAveraging','PCATemplate'};

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

fprintf('\nChecking optional functions:\n');
for i = 1:length(optionalMethods)
    fn = optionalMethods{i};
    fnPath = which(fn);
    if isempty(fnPath)
        fprintf('%s not found. It will be skipped if needed.\n',fn);
    else
        fprintf('%s -> %s\n',fn,fnPath);
    end
end

%% ========================================================================
% Method options
% ========================================================================

erpParams = struct();
erpParams.subtractWindow = artifactWindow;
erpParams.doBaselineCorrection = false;

ttaParams = struct();
ttaParams.window = artifactWindow;
ttaParams.computeFFT = false;

pcaParamsK10 = struct();
pcaParamsK10.artifactWindow = artifactWindow;
pcaParamsK10.numComponents = 10;
pcaParamsK10.removeMeanTemplate = true;
pcaParamsK10.taperEdgeMS = 2;
pcaParamsK10.doBaselineCorrection = false;

baseSMARTAOpts = struct();
baseSMARTAOpts.window = artifactWindow;
baseSMARTAOpts.computeFFT = false;
baseSMARTAOpts.pulseTimes = pulseTimes;
baseSMARTAOpts.stiFreq = smartaStimFreq;
baseSMARTAOpts.prePulse = smartaPrePulse;
baseSMARTAOpts.excludeSameTrial = true;

%% ========================================================================
% Search for best SMARTA/SMARTAFull representative
% ========================================================================

rows = struct([]);
rowCounter = 0;

best = struct();
best.objective = Inf;

bestLiteByElectrode = struct();
bestFullByElectrode = struct();

for iElec = 1:length(electrodesToSearch)

    elecNum = electrodesToSearch(iElec);
    elecFile = fullfile(lfpFolder,['elec' num2str(elecNum) '.mat']);

    fprintf('\n============================================================\n');
    fprintf('Searching elec%d\n',elecNum);
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

    fftNoStim = computeLogFFT(dataNoStim,fftIdx);
    fftRaw = computeLogFFT(dataHighStim,fftIdx);

    rawDiag = computeDiagnostics(fftRaw,fftNoStim,freqMask,harmonicMask);

    %% ERP

    erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);
    erpCleaned = getCleanedData(erpOut);
    fftERP = computeLogFFT(erpCleaned,fftIdx);
    erpDiag = computeDiagnostics(fftERP,fftNoStim,freqMask,harmonicMask);

    %% TTA, optional

    hasTTA = ~isempty(which('TrialWiseTemplateAveraging'));
    fftTTA = [];
    ttaDiag = [];

    if hasTTA
        try
            ttaOut = TrialWiseTemplateAveraging(dataHighStim,timeVals,ttaParams);
            ttaCleaned = getCleanedData(ttaOut);
            fftTTA = computeLogFFT(ttaCleaned,fftIdx);
            ttaDiag = computeDiagnostics(fftTTA,fftNoStim,freqMask,harmonicMask);
        catch ME
            warning('TTA failed on elec%d: %s',elecNum,ME.message);
            hasTTA = false;
        end
    end

    %% PCA, optional

    hasPCA = ~isempty(which('PCATemplate'));
    fftPCA = [];
    pcaDiag = [];

    if hasPCA
        try
            pcaOut = PCATemplate(dataHighStim,timeVals,pcaParamsK10);
            pcaCleaned = getCleanedData(pcaOut);
            fftPCA = computeLogFFT(pcaCleaned,fftIdx);
            pcaDiag = computeDiagnostics(fftPCA,fftNoStim,freqMask,harmonicMask);
        catch ME
            warning('PCA failed on elec%d: %s',elecNum,ME.message);
            hasPCA = false;
        end
    end

    %% SMARTALite search

    bestLite.objective = Inf;

    for iK = 1:length(liteKList)

        K = liteKList(iK);

        opts = baseSMARTAOpts;
        opts.K = K;

        try
            out = SMARTALite(dataHighStim,timeVals,opts);
            cleaned = getCleanedData(out);
            fftClean = computeLogFFT(cleaned,fftIdx);

            diagOut = computeDiagnostics(fftClean,fftNoStim,freqMask,harmonicMask);
            objective = diagOut.totalError + harmonicWeight * diagOut.harmonicAboveError;

            rowCounter = rowCounter + 1;
            rows(rowCounter).electrode = elecNum;
            rows(rowCounter).method = string('SMARTALite');
            rows(rowCounter).K = K;
            rows(rowCounter).hpCutoff = NaN;
            rows(rowCounter).totalError = diagOut.totalError;
            rows(rowCounter).harmonicAboveError = diagOut.harmonicAboveError;
            rows(rowCounter).harmonicSuppressionPct = 100 * ...
                (rawDiag.harmonicAboveError - diagOut.harmonicAboveError) / rawDiag.harmonicAboveError;
            rows(rowCounter).objective = objective;

            fprintf('SMARTALite K=%d: total %.4f, harmonic %.4f, objective %.4f\n', ...
                K,diagOut.totalError,diagOut.harmonicAboveError,objective);

            if objective < bestLite.objective
                bestLite.method = 'SMARTALite';
                bestLite.K = K;
                bestLite.hpCutoff = NaN;
                bestLite.fft = fftClean;
                bestLite.diag = diagOut;
                bestLite.objective = objective;
                bestLite.cleaned = cleaned;
            end

            if objective < best.objective
                best.elecNum = elecNum;
                best.method = 'SMARTALite';
                best.K = K;
                best.hpCutoff = NaN;
                best.fft = fftClean;
                best.diag = diagOut;
                best.objective = objective;

                best.dataNoStim = dataNoStim;
                best.dataHighStim = dataHighStim;
                best.fftNoStim = fftNoStim;
                best.fftRaw = fftRaw;
                best.fftERP = fftERP;
                best.fftTTA = fftTTA;
                best.fftPCA = fftPCA;
                best.rawDiag = rawDiag;
                best.erpDiag = erpDiag;
                best.ttaDiag = ttaDiag;
                best.pcaDiag = pcaDiag;
                best.hasTTA = hasTTA;
                best.hasPCA = hasPCA;
            end

        catch ME
            warning('SMARTALite K=%d failed on elec%d: %s',K,elecNum,ME.message);
        end
    end

    bestLiteByElectrode(iElec).electrode = elecNum;
    bestLiteByElectrode(iElec).result = bestLite;

    %% SMARTAFull search

    bestFull.objective = Inf;

    for iK = 1:length(fullKList)

        K = fullKList(iK);

        for iHp = 1:length(fullHpCutoffs)

            hpCutoff = fullHpCutoffs(iHp);

            opts = baseSMARTAOpts;
            opts.K = K;
            opts.hpCutoff = hpCutoff;
            opts.shrinkLoss = 'fro';
            opts.shrinkKL = 10;
            opts.shrinkKH = 15;

            try
                out = SMARTAFull(dataHighStim,timeVals,opts);
                cleaned = getCleanedData(out);
                fftClean = computeLogFFT(cleaned,fftIdx);

                diagOut = computeDiagnostics(fftClean,fftNoStim,freqMask,harmonicMask);
                objective = diagOut.totalError + harmonicWeight * diagOut.harmonicAboveError;

                rowCounter = rowCounter + 1;
                rows(rowCounter).electrode = elecNum;
                rows(rowCounter).method = string('SMARTAFull');
                rows(rowCounter).K = K;
                rows(rowCounter).hpCutoff = hpCutoff;
                rows(rowCounter).totalError = diagOut.totalError;
                rows(rowCounter).harmonicAboveError = diagOut.harmonicAboveError;
                rows(rowCounter).harmonicSuppressionPct = 100 * ...
                    (rawDiag.harmonicAboveError - diagOut.harmonicAboveError) / rawDiag.harmonicAboveError;
                rows(rowCounter).objective = objective;

                fprintf('SMARTAFull K=%d hp=%d: total %.4f, harmonic %.4f, objective %.4f\n', ...
                    K,hpCutoff,diagOut.totalError,diagOut.harmonicAboveError,objective);

                if objective < bestFull.objective
                    bestFull.method = 'SMARTAFull';
                    bestFull.K = K;
                    bestFull.hpCutoff = hpCutoff;
                    bestFull.fft = fftClean;
                    bestFull.diag = diagOut;
                    bestFull.objective = objective;
                    bestFull.cleaned = cleaned;
                end

                if objective < best.objective
                    best.elecNum = elecNum;
                    best.method = 'SMARTAFull';
                    best.K = K;
                    best.hpCutoff = hpCutoff;
                    best.fft = fftClean;
                    best.diag = diagOut;
                    best.objective = objective;

                    best.dataNoStim = dataNoStim;
                    best.dataHighStim = dataHighStim;
                    best.fftNoStim = fftNoStim;
                    best.fftRaw = fftRaw;
                    best.fftERP = fftERP;
                    best.fftTTA = fftTTA;
                    best.fftPCA = fftPCA;
                    best.rawDiag = rawDiag;
                    best.erpDiag = erpDiag;
                    best.ttaDiag = ttaDiag;
                    best.pcaDiag = pcaDiag;
                    best.hasTTA = hasTTA;
                    best.hasPCA = hasPCA;
                end

            catch ME
                warning('SMARTAFull K=%d hp=%d failed on elec%d: %s', ...
                    K,hpCutoff,elecNum,ME.message);
            end
        end
    end

    bestFullByElectrode(iElec).electrode = elecNum;
    bestFullByElectrode(iElec).result = bestFull;
end

if isempty(fieldnames(best))
    error('No SMARTA result was generated.');
end

searchTable = struct2table(rows);

%% ========================================================================
% Print selected best result
% ========================================================================

fprintf('\n============================================================\n');
fprintf('BEST REPRESENTATIVE SMARTAFull RESULT\n');
fprintf('============================================================\n');
fprintf('Electrode              : elec%d\n',best.elecNum);
fprintf('Best method            : %s\n',best.method);
fprintf('Best K                 : %d\n',best.K);

if strcmp(best.method,'SMARTAFull')
    fprintf('Best hpCutoff          : %.0f Hz\n',best.hpCutoff);
end

fprintf('Raw total FFT error    : %.4f\n',best.rawDiag.totalError);
fprintf('ERP total FFT error    : %.4f\n',best.erpDiag.totalError);

if best.hasTTA
    fprintf('TTA total FFT error    : %.4f\n',best.ttaDiag.totalError);
end

if best.hasPCA
    fprintf('PCA K10 total FFT error: %.4f\n',best.pcaDiag.totalError);
end

fprintf('SMARTA total FFT error : %.4f\n',best.diag.totalError);
fprintf('Raw harmonic error     : %.4f\n',best.rawDiag.harmonicAboveError);
fprintf('SMARTA harmonic error  : %.4f\n',best.diag.harmonicAboveError);
fprintf('Harmonic suppression   : %.2f %%\n', ...
    100*(best.rawDiag.harmonicAboveError - best.diag.harmonicAboveError) / best.rawDiag.harmonicAboveError);

%% ========================================================================
% Plot 1: Best SMARTA-focused FFT figure
% ========================================================================

fig1 = figure('Name','Best SMARTA FFT representative example', ...
    'Color','w','Position',[100 100 1050 600]);

plot(freqAxis(plotMask),best.fftNoStim(plotMask),'k','LineWidth',2.2);
hold on;
plot(freqAxis(plotMask),best.fftRaw(plotMask),'Color',[0.75 0.75 0.75],'LineWidth',1.6);
plot(freqAxis(plotMask),best.fftERP(plotMask),'LineWidth',1.5);

legendEntries = {'No-stim reference','Raw high-stim','ERPSubtraction'};

if best.hasTTA
    plot(freqAxis(plotMask),best.fftTTA(plotMask),'LineWidth',1.5);
    legendEntries{end+1} = 'TTA';
end

plot(freqAxis(plotMask),best.fft(plotMask),'LineWidth',2.4);

if strcmp(best.method,'SMARTAFull')
    legendEntries{end+1} = sprintf('SMARTAFull K=%d, hp=%d Hz',best.K,best.hpCutoff);
else
    legendEntries{end+1} = sprintf('SMARTALite K=%d',best.K);
end

xlim(freqRangeForPlot);
xlabel('Frequency (Hz)');
ylabel('log_{10} mean FFT magnitude (a.u.)');
title(sprintf('Representative SMARTAFull FFT result: elec%d',best.elecNum));
legend(legendEntries,'Location','best');
grid on;
box off;

if saveFigures
    saveAllFigureFormats(fig1,figFolder, ...
        sprintf('demo27_elec%d_best_SMARTA_fft',best.elecNum));
end

%% ========================================================================
% Plot 2: SMARTA vs PCA optional comparison
% ========================================================================

fig2 = figure('Name','Best SMARTA compared with baselines', ...
    'Color','w','Position',[100 100 1050 600]);

plot(freqAxis(plotMask),best.fftNoStim(plotMask),'k','LineWidth',2.2);
hold on;
plot(freqAxis(plotMask),best.fftRaw(plotMask),'Color',[0.75 0.75 0.75],'LineWidth',1.6);
plot(freqAxis(plotMask),best.fftERP(plotMask),'LineWidth',1.5);

legendEntries = {'No-stim reference','Raw high-stim','ERPSubtraction'};

if best.hasTTA
    plot(freqAxis(plotMask),best.fftTTA(plotMask),'LineWidth',1.5);
    legendEntries{end+1} = 'TTA';
end

if showPCAInFinalFigure && best.hasPCA
    plot(freqAxis(plotMask),best.fftPCA(plotMask),'LineWidth',1.5);
    legendEntries{end+1} = 'PCATemplate K=10';
end

plot(freqAxis(plotMask),best.fft(plotMask),'LineWidth',2.4);

if strcmp(best.method,'SMARTAFull')
    legendEntries{end+1} = sprintf('SMARTAFull K=%d, hp=%d Hz',best.K,best.hpCutoff);
else
    legendEntries{end+1} = sprintf('SMARTALite K=%d',best.K);
end

xlim(freqRangeForPlot);
xlabel('Frequency (Hz)');
ylabel('log_{10} mean FFT magnitude (a.u.)');
title(sprintf('SMARTA-focused FFT comparison: elec%d',best.elecNum));
legend(legendEntries,'Location','best');
grid on;
box off;

if saveFigures
    saveAllFigureFormats(fig2,figFolder, ...
        sprintf('demo27_elec%d_SMARTA_baseline_comparison',best.elecNum));
end

%% ========================================================================
% Plot 3: all searched SMARTA results for selected best electrode
% ========================================================================

idxBestElec = searchTable.electrode == best.elecNum;
TbestElec = sortrows(searchTable(idxBestElec,:),{'method','objective'});

fprintf('\nAll SMARTA search results for selected electrode:\n');
disp(TbestElec);

if saveFigures
    writetable(searchTable,fullfile(figFolder,'demo27_SMARTA_search_all_results.csv'));
    writetable(TbestElec,fullfile(figFolder, ...
        sprintf('demo27_elec%d_SMARTA_search_results.csv',best.elecNum)));
end

fprintf('\nSaved demo_27 figures and tables in:\n%s\n',figFolder);

%% ========================================================================
% Local helper functions
% ========================================================================

function y = computeLogFFT(data,fftIdx)

    x = data(:,fftIdx);
    fftVals = abs(fft(x,[],2));
    y = log10(mean(fftVals,1) + eps);
end

function diagOut = computeDiagnostics(fftMethod,fftReference,freqMask,harmonicMask)

    diffVals = fftMethod - fftReference;

    diffMetric = diffVals(freqMask);
    harmonicDiff = diffVals(harmonicMask);

    aboveVals = max(harmonicDiff,0);
    belowVals = max(-harmonicDiff,0);

    diagOut = struct();
    diagOut.totalError = norm(diffMetric);
    diagOut.harmonicAboveError = norm(aboveVals);
    diagOut.harmonicBelowError = norm(belowVals);
    diagOut.meanSignedDiff = mean(diffMetric);
end

function harmonicMask = buildHarmonicMask(freqAxis,harmonicsToUse,halfBandwidthHz,freqRangeForMetric)

    harmonicMask = false(size(freqAxis));

    for iH = 1:length(harmonicsToUse)
        h = harmonicsToUse(iH);

        if h < freqRangeForMetric(1) || h > freqRangeForMetric(2)
            continue;
        end

        harmonicMask = harmonicMask | abs(freqAxis - h) <= halfBandwidthHz;
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