clear; close all; clc;

%% demo_23_FinalPlotsAndTables
%
% Purpose:
% This script loads the latest final validation metric files from demo_22,
% generates summary tables, creates validation figures, and exports
% electrode-wise comparison results.
%
% Expected input files:
%   demo22_final_validation_summary_*.csv
%   demo22_final_validation_stats_*.csv
%   demo22_final_validation_best_counts_*.csv
%   demo22_final_validation_long_*.csv
%
% Main selected method in this script:
%   SMARTALite_Ensemble_K3K5
%
% Main comparison methods:
%   ERPSubtraction
%   PCATemplate_K10
%   SMARTALite_K3
%   SMARTALite_K5

%% ========================================================================
% Settings
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

metricsFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','metrics');

figFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','figures','demo23_FinalPlotsAndTables');

tableFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal','results','tables','demo23_FinalPlotsAndTables');

if ~exist(figFolder,'dir')
    mkdir(figFolder);
end

if ~exist(tableFolder,'dir')
    mkdir(tableFolder);
end

methodOrder = { ...
    'RawHighStim', ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_Ensemble_K3K5'};

methodLabels = { ...
    'Raw high-stim', ...
    'ERP subtraction', ...
    'PCA template K=10', ...
    'SMARTALite K=3', ...
    'SMARTALite K=5', ...
    'SMARTALite ensemble K3K5'};

mainMethods = { ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_Ensemble_K3K5'};

mainLabels = { ...
    'ERP subtraction', ...
    'PCA template K=10', ...
    'SMARTALite K=3', ...
    'SMARTALite K=5', ...
    'SMARTALite ensemble K3K5'};

timestampString = datestr(now,'yyyymmdd_HHMMSS');

%% ========================================================================
% Load latest demo_22 files
% ========================================================================

summaryFile = getLatestFile(metricsFolder,'demo22_final_validation_summary_*.csv');
statsFile   = getLatestFile(metricsFolder,'demo22_final_validation_stats_*.csv');
bestFile    = getLatestFile(metricsFolder,'demo22_final_validation_best_counts_*.csv');
longFile    = getLatestFile(metricsFolder,'demo22_final_validation_long_*.csv');

fprintf('\nUsing demo_22 files:\n');
fprintf('Summary: %s\n',summaryFile);
fprintf('Stats  : %s\n',statsFile);
fprintf('Best   : %s\n',bestFile);
fprintf('Long   : %s\n',longFile);

summaryTable = readtable(summaryFile);
statsTable = readtable(statsFile);
bestTable = readtable(bestFile);
longTable = readtable(longFile);

summaryTable.method = string(summaryTable.method);
statsTable.mainMethod = string(statsTable.mainMethod);
statsTable.comparisonMethod = string(statsTable.comparisonMethod);
statsTable.metric = string(statsTable.metric);

%% ========================================================================
% Reorder summary table
% ========================================================================

orderedSummary = table();

for i = 1:length(methodOrder)
    idx = summaryTable.method == string(methodOrder{i});
    if any(idx)
        orderedSummary = [orderedSummary; summaryTable(idx,:)];
    end
end

mainSummary = table();

for i = 1:length(mainMethods)
    idx = orderedSummary.method == string(mainMethods{i});
    if any(idx)
        mainSummary = [mainSummary; orderedSummary(idx,:)];
    end
end

orderedLabels = getLabelsForMethods(orderedSummary.method, methodOrder, methodLabels);
mainLabelsPresent = getLabelsForMethods(mainSummary.method, mainMethods, mainLabels);

fprintf('\nOrdered final summary table:\n');
disp(orderedSummary);

%% ========================================================================
% Export summary table
% ========================================================================

summaryTableOut = table();

summaryTableOut.Method = orderedLabels(:);
summaryTableOut.MeanTotalFFTError = orderedSummary.meanTotalFFTError;
summaryTableOut.StdTotalFFTError = orderedSummary.stdTotalFFTError;
summaryTableOut.MeanBalancedScore = orderedSummary.meanBalancedScore;
summaryTableOut.StdBalancedScore = orderedSummary.stdBalancedScore;
summaryTableOut.MeanHarmonicExcessError = orderedSummary.meanHarmonicAboveError;
summaryTableOut.MeanHarmonicSuppressionPct = orderedSummary.meanHarmonicSuppressionPct;
summaryTableOut.MeanBelowReferenceError = orderedSummary.meanBelowReferenceError;
summaryTableOut.MeanFractionBelowReference = orderedSummary.meanFractionBelowReference;
summaryTableOut.MeanNoStimDistortion = orderedSummary.meanNoStimDistortion;

summaryOutFile = fullfile(tableFolder, ...
    ['demo23_summary_table_' timestampString '.csv']);

writetable(summaryTableOut,summaryOutFile);

fprintf('\nSaved summary table:\n%s\n',summaryOutFile);
disp(summaryTableOut);

%% ========================================================================
% Figure 1: Main spectral metrics
% ========================================================================

fig1 = figure('Name','Final spectral validation metrics','Color','w');

barData = [ ...
    mainSummary.meanTotalFFTError, ...
    mainSummary.meanBalancedScore, ...
    mainSummary.meanHarmonicAboveError];

bar(barData);

xticks(1:height(mainSummary));
xticklabels(mainLabelsPresent);
xtickangle(25);

ylabel('Metric value');
title('Final validation: spectral artifact-removal metrics');
legend({'Total FFT error','Balanced score','Harmonic excess error'}, ...
    'Location','northwest');

grid on;
set(gca,'FontSize',11);

saveas(fig1,fullfile(figFolder,'demo23_final_spectral_metrics.png'));
savefig(fig1,fullfile(figFolder,'demo23_final_spectral_metrics.fig'));

%% ========================================================================
% Figure 2: Total FFT error with standard deviation
% ========================================================================

fig2 = figure('Name','Total FFT error','Color','w');

bar(mainSummary.meanTotalFFTError);
hold on;

errorbar(1:height(mainSummary), ...
    mainSummary.meanTotalFFTError, ...
    mainSummary.stdTotalFFTError, ...
    'k.','LineWidth',1.2);

xticks(1:height(mainSummary));
xticklabels(mainLabelsPresent);
xtickangle(25);

ylabel('Total FFT error vs no-stim reference');
title('Mean total FFT error across good V1 electrodes');
grid on;
set(gca,'FontSize',11);

saveas(fig2,fullfile(figFolder,'demo23_total_fft_error.png'));
savefig(fig2,fullfile(figFolder,'demo23_total_fft_error.fig'));

%% ========================================================================
% Figure 3: Balanced score
% ========================================================================

fig3 = figure('Name','Balanced score','Color','w');

bar(mainSummary.meanBalancedScore);
hold on;

errorbar(1:height(mainSummary), ...
    mainSummary.meanBalancedScore, ...
    mainSummary.stdBalancedScore, ...
    'k.','LineWidth',1.2);

xticks(1:height(mainSummary));
xticklabels(mainLabelsPresent);
xtickangle(25);

ylabel('Balanced score');
title('Balanced spectral score including below-reference penalty');
grid on;
set(gca,'FontSize',11);

saveas(fig3,fullfile(figFolder,'demo23_balanced_score.png'));
savefig(fig3,fullfile(figFolder,'demo23_balanced_score.fig'));

%% ========================================================================
% Figure 4: Above-reference vs below-reference error
% ========================================================================

fig4 = figure('Name','Above vs below reference error','Color','w');

barData = [ ...
    mainSummary.meanAboveReferenceError, ...
    mainSummary.meanBelowReferenceError];

bar(barData);

xticks(1:height(mainSummary));
xticklabels(mainLabelsPresent);
xtickangle(25);

ylabel('Spectral error norm');
title('Residual excess power vs possible over-cleaning');
legend({'Above no-stim reference','Below no-stim reference'}, ...
    'Location','northwest');

grid on;
set(gca,'FontSize',11);

saveas(fig4,fullfile(figFolder,'demo23_above_below_reference_error.png'));
savefig(fig4,fullfile(figFolder,'demo23_above_below_reference_error.fig'));

%% ========================================================================
% Figure 5: Harmonic suppression
% ========================================================================

fig5 = figure('Name','Harmonic suppression','Color','w');

bar(mainSummary.meanHarmonicSuppressionPct);
hold on;

errorbar(1:height(mainSummary), ...
    mainSummary.meanHarmonicSuppressionPct, ...
    mainSummary.stdHarmonicSuppressionPct, ...
    'k.','LineWidth',1.2);

xticks(1:height(mainSummary));
xticklabels(mainLabelsPresent);
xtickangle(25);

ylabel('Harmonic suppression relative to raw (%)');
title('Suppression of 20 Hz stimulation harmonics');
grid on;
set(gca,'FontSize',11);

saveas(fig5,fullfile(figFolder,'demo23_harmonic_suppression.png'));
savefig(fig5,fullfile(figFolder,'demo23_harmonic_suppression.fig'));

%% ========================================================================
% Figure 6: No-stim distortion control
% ========================================================================

fig6 = figure('Name','No-stim distortion control','Color','w');

bar(mainSummary.meanNoStimDistortion);
hold on;

errorbar(1:height(mainSummary), ...
    mainSummary.meanNoStimDistortion, ...
    mainSummary.stdNoStimDistortion, ...
    'k.','LineWidth',1.2);

xticks(1:height(mainSummary));
xticklabels(mainLabelsPresent);
xtickangle(25);

ylabel('No-stim distortion');
title('Control: distortion when methods are applied to no-stim trials');
grid on;
set(gca,'FontSize',11);

saveas(fig6,fullfile(figFolder,'demo23_no_stim_distortion.png'));
savefig(fig6,fullfile(figFolder,'demo23_no_stim_distortion.fig'));

%% ========================================================================
% Figure 7: Best-method counts
% ========================================================================

candidateMethods = { ...
    'ERPSubtraction', ...
    'PCATemplate_K10', ...
    'SMARTALite_K3', ...
    'SMARTALite_K5', ...
    'SMARTALite_Ensemble_K3K5'};

candidateLabels = mainLabels;

countTotal = countBest(bestTable.bestTotalFFT,candidateMethods);
countBalanced = countBest(bestTable.bestBalanced,candidateMethods);
countHarmonic = countBest(bestTable.bestHarmonicSuppression,candidateMethods);

countData = [countTotal(:), countBalanced(:), countHarmonic(:)];

fig7 = figure('Name','Best method counts','Color','w');

bar(countData);

xticks(1:length(candidateMethods));
xticklabels(candidateLabels);
xtickangle(25);

ylabel('Number of electrodes');
title('Best-method counts across 31 good V1 electrodes');
legend({'Total FFT error','Balanced score','Harmonic excess error'}, ...
    'Location','northwest');

grid on;
set(gca,'FontSize',11);

saveas(fig7,fullfile(figFolder,'demo23_best_method_counts.png'));
savefig(fig7,fullfile(figFolder,'demo23_best_method_counts.fig'));

bestCountsTable = table( ...
    candidateLabels(:), ...
    countTotal(:), ...
    countBalanced(:), ...
    countHarmonic(:), ...
    'VariableNames',{'Method','BestTotalFFTCount','BestBalancedCount','BestHarmonicCount'} ...
);

bestCountsOutFile = fullfile(tableFolder, ...
    ['demo23_best_method_counts_' timestampString '.csv']);

writetable(bestCountsTable,bestCountsOutFile);

fprintf('\nSaved best-method counts table:\n%s\n',bestCountsOutFile);
disp(bestCountsTable);

%% ========================================================================
% Figure 8: Electrode-wise total FFT error
% ========================================================================

fig8 = figure('Name','Electrode-wise total FFT error','Color','w');

hold on;

plotMetricByElectrode(longTable,'PCATemplate_K10','totalFFTError','o-');
plotMetricByElectrode(longTable,'SMARTALite_K3','totalFFTError','o-');
plotMetricByElectrode(longTable,'SMARTALite_Ensemble_K3K5','totalFFTError','o-');

xlabel('Electrode');
ylabel('Total FFT error');
title('Electrode-wise total FFT error');
legend({'PCA template K=10','SMARTALite K=3','SMARTALite ensemble K3K5'}, ...
    'Location','best');

grid on;
set(gca,'FontSize',11);

saveas(fig8,fullfile(figFolder,'demo23_electrodewise_total_fft_error.png'));
savefig(fig8,fullfile(figFolder,'demo23_electrodewise_total_fft_error.fig'));

%% ========================================================================
% Export statistical comparison table
% ========================================================================

statsKeep = statsTable( ...
    statsTable.mainMethod == "SMARTALite_Ensemble_K3K5" & ...
    ismember(statsTable.metric,["totalFFTError","balancedScore","harmonicAboveError","trialSTDErrorVsNoStim"]),:);

statsOutFile = fullfile(tableFolder, ...
    ['demo23_ensemble_pairwise_statistics_' timestampString '.csv']);

writetable(statsKeep,statsOutFile);

fprintf('\nSaved pairwise statistics table:\n%s\n',statsOutFile);
disp(statsKeep);

%% ========================================================================
% Final printed summary
% ========================================================================

ensembleRow = orderedSummary(orderedSummary.method == "SMARTALite_Ensemble_K3K5",:);
k3Row = orderedSummary(orderedSummary.method == "SMARTALite_K3",:);
pcaRow = orderedSummary(orderedSummary.method == "PCATemplate_K10",:);

fprintf('\n============================================================\n');
fprintf('FINAL SUMMARY NUMBERS\n');
fprintf('============================================================\n');

fprintf('SMARTALite Ensemble K3K5:\n');
fprintf('  Mean total FFT error      = %.4f\n',ensembleRow.meanTotalFFTError);
fprintf('  Mean balanced score       = %.4f\n',ensembleRow.meanBalancedScore);
fprintf('  Mean harmonic excess      = %.4f\n',ensembleRow.meanHarmonicAboveError);
fprintf('  Harmonic suppression      = %.2f %%\n',ensembleRow.meanHarmonicSuppressionPct);
fprintf('  Mean below-reference err  = %.4f\n',ensembleRow.meanBelowReferenceError);
fprintf('  Mean no-stim distortion   = %.4f\n',ensembleRow.meanNoStimDistortion);

fprintf('\nSMARTALite K3:\n');
fprintf('  Mean total FFT error      = %.4f\n',k3Row.meanTotalFFTError);
fprintf('  Mean balanced score       = %.4f\n',k3Row.meanBalancedScore);
fprintf('  Mean harmonic excess      = %.4f\n',k3Row.meanHarmonicAboveError);

fprintf('\nPCATemplate K10:\n');
fprintf('  Mean total FFT error      = %.4f\n',pcaRow.meanTotalFFTError);
fprintf('  Mean balanced score       = %.4f\n',pcaRow.meanBalancedScore);
fprintf('  Mean harmonic excess      = %.4f\n',pcaRow.meanHarmonicAboveError);

fprintf('\nSaved final figures in:\n%s\n',figFolder);
fprintf('Saved final tables in:\n%s\n',tableFolder);

%% ========================================================================
% Local helper functions
% ========================================================================

function latestFile = getLatestFile(folderPath,pattern)

    files = dir(fullfile(folderPath,pattern));

    if isempty(files)
        error('No files found for pattern: %s',fullfile(folderPath,pattern));
    end

    [~,idx] = max([files.datenum]);
    latestFile = fullfile(files(idx).folder,files(idx).name);
end

function counts = countBest(bestColumn,candidateMethods)

    bestColumn = string(bestColumn);
    counts = zeros(length(candidateMethods),1);

    for i = 1:length(candidateMethods)
        counts(i) = sum(bestColumn == string(candidateMethods{i}));
    end
end

function labels = getLabelsForMethods(methodNames, methodOrder, methodLabels)

    methodNames = string(methodNames);
    labels = strings(length(methodNames),1);

    for i = 1:length(methodNames)
        idx = strcmp(methodOrder, char(methodNames(i)));

        if any(idx)
            labels(i) = string(methodLabels{idx});
        else
            labels(i) = methodNames(i);
        end
    end

    labels = cellstr(labels);
end

function plotMetricByElectrode(longTable,methodName,metricName,lineSpec)

    methodName = string(methodName);

    idx = string(longTable.method) == methodName;
    T = longTable(idx,:);

    if isempty(T)
        warning('No rows found for method: %s',methodName);
        return;
    end

    [~,sortIdx] = sort(T.electrode);
    T = T(sortIdx,:);

    plot(T.electrode,T.(metricName),lineSpec,'LineWidth',1.3);
end