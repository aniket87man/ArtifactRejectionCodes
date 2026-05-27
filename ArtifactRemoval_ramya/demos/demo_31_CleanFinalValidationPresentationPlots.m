%% demo_31_CleanFinalValidationPresentationPlots
%
% Clean plots from demo_29 final validation results.
%


clear; close all; clc;

%% ========================================================================
% User paths
% ========================================================================

folderSourceString = 'C:\Users\RAMYA\Documents\MATLAB\programs';

metricsFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal', ...
    'results', ...
    'metrics', ...
    'demo29_FinalValidation_WithSMARTAFull');

figFolder = fullfile(folderSourceString, ...
    'icms_artifact_removal', ...
    'results', ...
    'figures', ...
    'demo31_CleanPresentationPlots');

if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

summaryFile = fullfile(metricsFolder, 'demo29_final_validation_summary.csv');

if ~exist(summaryFile, 'file')
    error('Summary file not found:\n%s\nRun demo_29_FinalValidation first.', summaryFile);
end

%% ========================================================================
% Load summary table
% ========================================================================

T = readtable(summaryFile);

fprintf('\nLoaded summary table:\n%s\n', summaryFile);
disp(T);

% Convert method names to string array
methodNames = string(T.method);

% Display names for plots
displayNames = methodNames;
displayNames(displayNames == "RawHighStim") = "Raw high-stim";
displayNames(displayNames == "ERPSubtraction") = "ERPSubtraction";
displayNames(displayNames == "PCATemplate_K10") = "PCATemplate K=10";
displayNames(displayNames == "SMARTALite_Ensemble_K3K5") = "SMARTALite ensemble K3K5";
displayNames(displayNames == "SMARTAFull_K3_hp100") = "SMARTAFull K=3, hp=100 Hz";

%% ========================================================================
% Plot settings
% ========================================================================

fontSize = 14;
titleFontSize = 16;
barWidth = 0.65;

%% ========================================================================
% Figure 1: Total FFT error
% Raw high-stim is included because this is a real contaminated baseline.
% ========================================================================

plotMask = true(size(methodNames));

makeBarPlot( ...
    T.meanTotalFFTError, ...
    getOptionalColumn(T, 'stdTotalFFTError'), ...
    displayNames, ...
    plotMask, ...
    'Final validation: total FFT error (lower is better)', ...
    'Mean total FFT error', ...
    fullfile(figFolder, 'demo31_mean_total_fft_error_clean'), ...
    fontSize, titleFontSize, barWidth);

%% ========================================================================
% Figure 2: Harmonic excess error
% Raw high-stim is included because this shows stimulation artifact strength.
% ========================================================================

plotMask = true(size(methodNames));

makeBarPlot( ...
    T.meanHarmonicExcessError, ...
    getOptionalColumn(T, 'stdHarmonicExcessError'), ...
    displayNames, ...
    plotMask, ...
    'Final validation: stimulation-harmonic excess error (lower is better)', ...
    'Mean harmonic excess error', ...
    fullfile(figFolder, 'demo31_mean_harmonic_excess_error_clean'), ...
    fontSize, titleFontSize, barWidth);

%% ========================================================================
% Figure 3: Balanced score
% Remove Raw high-stim because its balanced score is NaN.
% ========================================================================

plotMask = methodNames ~= "RawHighStim";

makeBarPlot( ...
    T.meanBalancedScore, ...
    getOptionalColumn(T, 'stdBalancedScore'), ...
    displayNames, ...
    plotMask, ...
    'Final validation: normalized balanced score (lower is better)', ...
    'Mean normalized balanced score', ...
    fullfile(figFolder, 'demo31_mean_balanced_score_clean'), ...
    fontSize, titleFontSize, barWidth);

%% ========================================================================
% Figure 4: No-stimulation negative-control FFT distortion
% Remove Raw high-stim because raw/no-cleaning gives zero distortion by definition.
% ========================================================================

plotMask = methodNames ~= "RawHighStim";

makeBarPlot( ...
    T.meanNoStimFFTDistortion, ...
    [], ...
    displayNames, ...
    plotMask, ...
    'Negative-control distortion on no-stimulation trials (lower is better)', ...
    'Mean no-stimulation FFT distortion', ...
    fullfile(figFolder, 'demo31_mean_no_stim_fft_distortion_clean'), ...
    fontSize, titleFontSize, barWidth);

fprintf('\nClean presentation plots saved in:\n%s\n', figFolder);

%% ========================================================================
% Local helper functions
% ========================================================================

function makeBarPlot(values, stdValues, displayNames, plotMask, titleStr, yLabelStr, outBase, fontSize, titleFontSize, barWidth)

    values = values(:);
    displayNames = displayNames(:);
    plotMask = plotMask(:);

    y = values(plotMask);
    labels = displayNames(plotMask);

    if ~isempty(stdValues)
        stdValues = stdValues(:);
        err = stdValues(plotMask);
    else
        err = [];
    end

    fig = figure('Color', 'w', 'Position', [100 100 1200 650]);

    b = bar(y, barWidth);
    b.LineWidth = 1.0;

    hold on;

    if ~isempty(err)
        x = 1:numel(y);
        errorbar(x, y, err, 'k.', 'LineWidth', 1.2, 'CapSize', 12);
    end

    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(25);

    ylabel(yLabelStr, 'FontSize', fontSize);
    title(titleStr, 'FontSize', titleFontSize, 'FontWeight', 'bold');

    set(gca, 'FontSize', fontSize);
    grid on;
    box off;

    ymax = max(y + max([zeros(size(y)), getErrForLimit(err, y)], [], 2), [], 'omitnan');

    if isempty(ymax) || isnan(ymax) || ymax <= 0
        ymax = 1;
    end

    ylim([0 1.18*ymax]);

    % Save
    pngFile = [outBase '.png'];
    pdfFile = [outBase '.pdf'];
    figFile = [outBase '.fig'];

    saveas(fig, pngFile);
    savefig(fig, figFile);

    try
        exportgraphics(fig, pdfFile, 'ContentType', 'vector');
    catch
        exportgraphics(fig, pdfFile, 'ContentType', 'image');
    end

    fprintf('\nSaved:\n%s\n%s\n%s\n', pngFile, pdfFile, figFile);
end

function errForLimit = getErrForLimit(err, y)
    if isempty(err)
        errForLimit = zeros(size(y));
    else
        errForLimit = err;
        errForLimit(isnan(errForLimit)) = 0;
    end
end

function col = getOptionalColumn(T, colName)
    if ismember(colName, T.Properties.VariableNames)
        col = T.(colName);
    else
        col = [];
    end
end