%% demo_13_PCATemplateCrossValidation
% Cross-validate PCATemplate for PCA-based artifact removal.
%
% This version:
% 1. Uses elec7 as the representative analyzed electrode.
% 2. Uses the same artifact/FFT window: 0--0.4 s.
% 3. Tests K = 0, 1, 2, 3, 5, 10.
% 4. Compares PCATemplate only with ERPSubtraction.
% 5. Saves:
%       demo13_pca_cv_improvement_elec7.pdf
%       demo13_pca_fft_simplified_elec7.pdf

clear; close all; clc;

%% ========================================================================
% User settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

subjectName   = 'dona';
gridType      = 'Microelectrode';
expDate       = '290825';
protocolName  = 'GRF_001';
electrodeName = 'elec7';

% Conditions
noStimCondition   = {1,1,1,5,5,4};
highStimCondition = {7,1,1,5,5,4};

% PCA settings
artifactWindow = [0 0.4];
kList = [0 1 2 3 5 10];

% Cross-validation settings
numFolds = 5;
numRepeats = 10;
rngSeed = 1;

% FFT settings
fftWindow = [0 0.4];
freqRangeForMetric = [0 200];

%% ========================================================================
% Paths and load data
% ========================================================================

baseFolder = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

lfpFile     = fullfile(baseFolder,'segmentedData','lfp',[electrodeName '.mat']);
lfpInfoFile = fullfile(baseFolder,'segmentedData','lfp','lfpInfo.mat');
paramFile   = fullfile(baseFolder,'extractedData','parameterCombinations.mat');

figFolder = fullfile(pwd,'figures');
if ~exist(figFolder,'dir')
    mkdir(figFolder);
end

D = load(lfpFile);
I = load(lfpInfoFile);
P = load(paramFile);

analogData = D.analogData;
timeVals = I.timeVals;
parameterCombinations = P.parameterCombinations;

%% ========================================================================
% Select trials
% ========================================================================

noStimTrials = parameterCombinations{noStimCondition{:}};
highStimTrials = parameterCombinations{highStimCondition{:}};

dataNoStim = analogData(noStimTrials,:);
dataHighStim = analogData(highStimTrials,:);

fprintf('No-stim trials   : %d\n',size(dataNoStim,1));
fprintf('High-stim trials : %d\n',size(dataHighStim,1));

numHighTrials = size(dataHighStim,1);

%% ========================================================================
% FFT indices
% ========================================================================

% This gives 800 samples: 0.0005 to 0.4000 s
fftIdx = find(timeVals > fftWindow(1) & timeVals <= fftWindow(2));

fprintf('\nFFT window has %d samples.\n',length(fftIdx));
fprintf('Actual FFT window: %.4f to %.4f s\n',timeVals(fftIdx(1)),timeVals(fftIdx(end)));
fprintf('Frequency resolution: %.2f Hz\n',1/median(diff(timeVals))/length(fftIdx));

%% ========================================================================
% No-stim reference FFT
% ========================================================================

fftNoStim = compute_fft_summary_local(dataNoStim,timeVals,fftIdx);

freqMask = fftNoStim.freqAxis >= freqRangeForMetric(1) & ...
           fftNoStim.freqAxis <= freqRangeForMetric(2);

%% ========================================================================
% Reference full-data method: ERPSubtraction
% ========================================================================

erpParams = struct();
erpParams.subtractWindow = artifactWindow;
erpParams.doBaselineCorrection = false;

erpOut = ERPSubtraction(dataHighStim,timeVals,erpParams);

fftHighRawFull = compute_fft_summary_local(dataHighStim,timeVals,fftIdx);
fftERPClean = compute_fft_summary_local(erpOut.cleanedData,timeVals,fftIdx);

rawFullError = norm(fftHighRawFull.logMeanMagnitude(freqMask) - ...
                    fftNoStim.logMeanMagnitude(freqMask));

erpError = norm(fftERPClean.logMeanMagnitude(freqMask) - ...
                fftNoStim.logMeanMagnitude(freqMask));

erpImprovement = 100 * (rawFullError - erpError) / rawFullError;

fprintf('\nReference full-data method results:\n');
fprintf('Raw high-stim FFT error vs no-stim = %.4f\n',rawFullError);
fprintf('ERPSubtraction error               = %.4f\n',erpError);
fprintf('ERPSubtraction improvement         = %.2f %%\n',erpImprovement);

%% ========================================================================
% Cross-validation
% ========================================================================

rng(rngSeed);

cvRows = struct([]);
rowCounter = 0;

fprintf('\nRunning PCATemplate cross-validation...\n');
fprintf('numFolds = %d, numRepeats = %d\n',numFolds,numRepeats);

for iRepeat = 1:numRepeats

    shuffledIdx = randperm(numHighTrials);
    foldIDs = zeros(numHighTrials,1);

    for i = 1:numHighTrials
        foldIDs(shuffledIdx(i)) = mod(i-1,numFolds) + 1;
    end

    for iFold = 1:numFolds

        testMask = foldIDs == iFold;
        trainMask = ~testMask;

        trainData = dataHighStim(trainMask,:);
        testData = dataHighStim(testMask,:);

        fftTestRaw = compute_fft_summary_local(testData,timeVals,fftIdx);

        rawTestError = norm(fftTestRaw.logMeanMagnitude(freqMask) - ...
                            fftNoStim.logMeanMagnitude(freqMask));

        for iK = 1:length(kList)

            K = kList(iK);

            trainParams = struct();
            trainParams.artifactWindow = artifactWindow;
            trainParams.numComponents = K;
            trainParams.removeMeanTemplate = true;
            trainParams.taperEdgeMS = 2;
            trainParams.doBaselineCorrection = false;

            model = trainPCATemplateModel(trainData,timeVals,trainParams);

            [cleanTestData,artifactModelTest] = applyPCATemplateModel(testData,timeVals,model);

            fftCleanTest = compute_fft_summary_local(cleanTestData,timeVals,fftIdx);

            cleanTestError = norm(fftCleanTest.logMeanMagnitude(freqMask) - ...
                                  fftNoStim.logMeanMagnitude(freqMask));

            improvement = 100 * (rawTestError - cleanTestError) / rawTestError;

            rowCounter = rowCounter + 1;

            cvRows(rowCounter).repeat = iRepeat;
            cvRows(rowCounter).fold = iFold;
            cvRows(rowCounter).K = K;
            cvRows(rowCounter).numTrainTrials = size(trainData,1);
            cvRows(rowCounter).numTestTrials = size(testData,1);
            cvRows(rowCounter).rawTestError = rawTestError;
            cvRows(rowCounter).cleanTestError = cleanTestError;
            cvRows(rowCounter).improvement = improvement;
            cvRows(rowCounter).meanArtifactModelRMS = sqrt(mean(artifactModelTest(:).^2));

            fprintf('Repeat %02d/%02d, fold %d/%d, K=%2d -> raw %.4f, clean %.4f, improvement %.2f %%\n', ...
                iRepeat,numRepeats,iFold,numFolds,K,rawTestError,cleanTestError,improvement);

        end
    end
end

cvTable = struct2table(cvRows);

fprintf('\nCross-validation table preview:\n');
disp(cvTable(1:min(10,height(cvTable)),:));

%% ========================================================================
% Summarize cross-validation results
% ========================================================================

summaryRows = struct([]);

for iK = 1:length(kList)

    K = kList(iK);
    rowsK = cvTable(cvTable.K == K,:);

    summaryRows(iK).K = K;

    summaryRows(iK).meanRawTestError = mean(rowsK.rawTestError);
    summaryRows(iK).stdRawTestError = std(rowsK.rawTestError);

    summaryRows(iK).meanCleanTestError = mean(rowsK.cleanTestError);
    summaryRows(iK).stdCleanTestError = std(rowsK.cleanTestError);

    summaryRows(iK).meanImprovement = mean(rowsK.improvement);
    summaryRows(iK).stdImprovement = std(rowsK.improvement);

    summaryRows(iK).meanArtifactModelRMS = mean(rowsK.meanArtifactModelRMS);
    summaryRows(iK).stdArtifactModelRMS = std(rowsK.meanArtifactModelRMS);

end

summaryTable = struct2table(summaryRows);

fprintf('\nPCATemplate cross-validation summary:\n');
disp(summaryTable);

[bestMeanError,bestIdx] = min(summaryTable.meanCleanTestError);
bestK = summaryTable.K(bestIdx);
bestMeanImprovement = summaryTable.meanImprovement(bestIdx);

fprintf('\nBest K by mean held-out FFT error:\n');
fprintf('Best K = %d\n',bestK);
fprintf('Mean held-out clean error = %.4f\n',bestMeanError);
fprintf('Mean held-out improvement = %.2f %%\n',bestMeanImprovement);

%% Save CV tables
cvFile = fullfile(figFolder,['demo13_pca_cv_table_' electrodeName '.csv']);
summaryFile = fullfile(figFolder,['demo13_pca_cv_summary_' electrodeName '.csv']);

writetable(cvTable,cvFile);
writetable(summaryTable,summaryFile);

fprintf('\nSaved CV table:\n%s\n',cvFile);
fprintf('Saved CV summary:\n%s\n',summaryFile);

%% ========================================================================
% Train final model on all high-stim trials using best cross-validated K
% ========================================================================

bestParams = struct();
bestParams.artifactWindow = artifactWindow;
bestParams.numComponents = bestK;
bestParams.removeMeanTemplate = true;
bestParams.taperEdgeMS = 2;
bestParams.doBaselineCorrection = false;

bestModelFull = trainPCATemplateModel(dataHighStim,timeVals,bestParams);
[bestPCACleanFull,bestPCAArtifactFull] = applyPCATemplateModel(dataHighStim,timeVals,bestModelFull);

fftBestPCAFull = compute_fft_summary_local(bestPCACleanFull,timeVals,fftIdx);

bestFullError = norm(fftBestPCAFull.logMeanMagnitude(freqMask) - ...
                     fftNoStim.logMeanMagnitude(freqMask));

bestFullImprovement = 100 * (rawFullError - bestFullError) / rawFullError;

fprintf('\nFinal full-data PCATemplate using best CV K:\n');
fprintf('Best CV K = %d\n',bestK);
fprintf('Full-data PCA error = %.4f\n',bestFullError);
fprintf('Full-data PCA improvement = %.2f %%\n',bestFullImprovement);

%% ========================================================================
% Figure 1: Held-out error vs K
% ========================================================================

hFigError = figure('Name','PCATemplate CV: Held-out error vs K');

errorbar(summaryTable.K, ...
         summaryTable.meanCleanTestError, ...
         summaryTable.stdCleanTestError, ...
         'o-','LineWidth',1.2);

hold on;
yline(erpError,'b--','LineWidth',1.2);

xlabel('Number of PCA components, K');
ylabel('Held-out FFT error vs no-stim');
title('PCATemplate component selection by cross-validation');
legend({'PCATemplate CV','ERPSubtraction'},'Location','best');
grid on;
box off;

outFileError = fullfile(figFolder,['demo13_pca_cv_error_' electrodeName '.pdf']);
exportgraphics(hFigError,outFileError,'ContentType','vector');
fprintf('Saved PCA CV error figure:\n%s\n',outFileError);

%% ========================================================================
% Figure 2: Held-out improvement vs K
% ========================================================================

hFigCV = figure('Name','PCATemplate CV: Held-out improvement vs K');

errorbar(summaryTable.K, ...
         summaryTable.meanImprovement, ...
         summaryTable.stdImprovement, ...
         'o-','LineWidth',1.2);

hold on;
yline(erpImprovement,'b--','LineWidth',1.2);

xlabel('Number of PCA components, K');
ylabel('Held-out FFT-error improvement (%)');
title('PCATemplate component selection by cross-validation');
legend({'PCATemplate CV','ERPSubtraction'},'Location','best');
grid on;
box off;

outFileCV = fullfile(figFolder,['demo13_pca_cv_improvement_' electrodeName '.pdf']);
exportgraphics(hFigCV,outFileCV,'ContentType','vector');
fprintf('Saved PCA CV improvement figure:\n%s\n',outFileCV);

%% ========================================================================
% Figure 3: Simplified FFT comparison using best CV K
% ========================================================================

hFigFFT = figure('Name','PCATemplate CV: Simplified FFT comparison');

plot(fftNoStim.freqAxis,fftNoStim.logMeanMagnitude,'k','LineWidth',1.5);
hold on;
plot(fftHighRawFull.freqAxis,fftHighRawFull.logMeanMagnitude,'r','LineWidth',1.3);
plot(fftERPClean.freqAxis,fftERPClean.logMeanMagnitude,'b','LineWidth',1.3);
plot(fftBestPCAFull.freqAxis,fftBestPCAFull.logMeanMagnitude,'c','LineWidth',1.5);

xlim([0 200]);
xlabel('Frequency (Hz)');
ylabel('log_{10} mean FFT magnitude (a.u.)');
title('Frequency-domain comparison of PCATemplate and ERPSubtraction');

legend({'No-stimulation reference', ...
        'Uncleaned 64 \muA', ...
        'ERPSubtraction', ...
        ['PCATemplate K = ' num2str(bestK)]}, ...
        'Location','best');

box off;

outFileFFT = fullfile(figFolder,['demo13_pca_fft_simplified_' electrodeName '.pdf']);
exportgraphics(hFigFFT,outFileFFT,'ContentType','vector');

fprintf('Saved simplified PCA FFT figure:\n%s\n',outFileFFT);

%% ========================================================================
% Final printed interpretation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('PCATemplate cross-validation summary\n');
fprintf('============================================================\n');

fprintf('Artifact window = [%.3f %.3f] s\n',artifactWindow(1),artifactWindow(2));
fprintf('FFT metric range = [%d %d] Hz\n',freqRangeForMetric(1),freqRangeForMetric(2));
fprintf('Number of folds = %d\n',numFolds);
fprintf('Number of repeats = %d\n',numRepeats);

fprintf('\nReference full-data method:\n');
fprintf('ERPSubtraction: error %.4f, improvement %.2f %%\n',erpError,erpImprovement);

fprintf('\nBest cross-validated PCATemplate:\n');
fprintf('Best K = %d\n',bestK);
fprintf('Mean held-out clean error = %.4f\n',bestMeanError);
fprintf('Mean held-out improvement = %.2f %%\n',bestMeanImprovement);

fprintf('\nFull-data result using best CV K:\n');
fprintf('Full-data PCA error = %.4f\n',bestFullError);
fprintf('Full-data PCA improvement = %.2f %%\n',bestFullImprovement);

fprintf('\nInterpretation:\n');
fprintf('Higher K values subtract more PCA components and are therefore more aggressive.\n');
fprintf('If a lower K performs similarly to the best K, the lower K should be preferred as the more conservative PCA baseline.\n');
fprintf('In this run, K=%d gave the lowest held-out FFT error among the tested values.\n',bestK);

fprintf('\ndemo_13_PCATemplateCrossValidation complete.\n');

%% ========================================================================
% Local helper functions
% ========================================================================

function fftOut = compute_fft_summary_local(dataIn,timeVals,fftIdx)

    X = dataIn(:,fftIdx);

    Fs = 1/median(diff(timeVals));
    N = length(fftIdx);

    freqAxis = (0:N-1) * (Fs/N);

    
    % log10(mean(abs(fft(trials)')))
    fftVals = fft(X');
    logMeanMagnitude = log10(mean(abs(fftVals)',1));

    fftOut = struct();
    fftOut.freqAxis = freqAxis;
    fftOut.logMeanMagnitude = logMeanMagnitude;

end

function model = trainPCATemplateModel(trainData,timeVals,params)

    if ~isfield(params,'artifactWindow') || isempty(params.artifactWindow)
        params.artifactWindow = [0 0.4];
    end

    if ~isfield(params,'numComponents') || isempty(params.numComponents)
        params.numComponents = 3;
    end

    if ~isfield(params,'removeMeanTemplate') || isempty(params.removeMeanTemplate)
        params.removeMeanTemplate = true;
    end

    if ~isfield(params,'taperEdgeMS') || isempty(params.taperEdgeMS)
        params.taperEdgeMS = 2;
    end

    artifactIdx = find(timeVals > params.artifactWindow(1) & ...
                       timeVals <= params.artifactWindow(2));

    if isempty(artifactIdx)
        error('No samples found in artifactWindow.');
    end

    Fs = 1/median(diff(timeVals));

    X = trainData(:,artifactIdx);

    meanTemplateWindow = mean(X,1);

    if params.removeMeanTemplate
        Xcentered = X - meanTemplateWindow;
    else
        Xcentered = X;
        meanTemplateWindow = zeros(size(meanTemplateWindow));
    end

    [~,~,V] = svd(Xcentered,'econ');

    maxComponents = min([params.numComponents,size(V,2)]);

    if maxComponents < 1
        pcTemplates = zeros(0,length(artifactIdx));
    else
        pcTemplates = V(:,1:maxComponents)';
    end

    edgeSamples = round((params.taperEdgeMS/1000) * Fs);

    model = struct();
    model.artifactIdx = artifactIdx;
    model.meanTemplateWindow = meanTemplateWindow;
    model.pcTemplates = pcTemplates;
    model.numComponents = maxComponents;
    model.params = params;
    model.edgeSamples = edgeSamples;

end

function [cleanData,artifactModel] = applyPCATemplateModel(testData,timeVals,model)

    numTrials = size(testData,1);
    numTime = size(testData,2);

    artifactIdx = model.artifactIdx;

    Xtest = testData(:,artifactIdx);

    K = model.numComponents;

    if K < 1
        artifactModelWindow = repmat(model.meanTemplateWindow,numTrials,1);
    else
        XcenteredTest = Xtest - model.meanTemplateWindow;

        pcTemplates = model.pcTemplates;
        scoresTest = XcenteredTest * pcTemplates';
        pcaModelCentered = scoresTest * pcTemplates;

        artifactModelWindow = repmat(model.meanTemplateWindow,numTrials,1) + pcaModelCentered;
    end

    artifactModelWindow = applyWindowTaperLocal(artifactModelWindow,model.edgeSamples);

    artifactModel = zeros(numTrials,numTime);
    artifactModel(:,artifactIdx) = artifactModelWindow;

    cleanData = testData - artifactModel;

end

function Y = applyWindowTaperLocal(X,edgeSamples)

    Y = X;

    if edgeSamples <= 0
        return;
    end

    n = size(X,2);
    edgeSamples = min(edgeSamples,floor(n/2));

    if edgeSamples <= 0
        return;
    end

    t = linspace(0,pi/2,edgeSamples);
    leftTaper = sin(t).^2;
    rightTaper = fliplr(leftTaper);

    Y(:,1:edgeSamples) = Y(:,1:edgeSamples) .* leftTaper;
    Y(:,end-edgeSamples+1:end) = Y(:,end-edgeSamples+1:end) .* rightTaper;

end