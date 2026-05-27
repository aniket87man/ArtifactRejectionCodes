function out = ERPScaled(trialData,timeVals,params)
% ERPScaled
%
% ERP subtraction with trial-wise amplitude scaling.
%
% The ERP template is computed across trials. For each trial, a scale factor
% alpha is estimated by least-squares fitting the ERP template to that trial
% inside scaleWindow. The scaled template is then subtracted only inside
% subtractWindow.
%
% Final window convention:
%   [0 0.4] means timeVals > 0 and timeVals <= 0.4
%   This gives 0.0005 to 0.4000 s, i.e., 800 samples at 2000 Hz.

if nargin < 3
    params = struct();
end

if ~isfield(params,'subtractWindow') || isempty(params.subtractWindow)
    params.subtractWindow = [timeVals(1)-eps timeVals(end)];
end

if ~isfield(params,'scaleWindow') || isempty(params.scaleWindow)
    params.scaleWindow = params.subtractWindow;
end

if ~isfield(params,'scaleBounds') || isempty(params.scaleBounds)
    params.scaleBounds = [0 2];
end

if ~isfield(params,'doBaselineCorrection') || isempty(params.doBaselineCorrection)
    params.doBaselineCorrection = false;
end

if ~isfield(params,'baselineWindow') || isempty(params.baselineWindow)
    params.baselineWindow = [-0.7 -0.2];
end

if size(trialData,2) ~= length(timeVals)
    error('trialData must be trials x time, and length(timeVals) must match size(trialData,2).');
end

dataIn = trialData;

%% Optional baseline correction
if params.doBaselineCorrection
    baselineIdx = timeVals > params.baselineWindow(1) & ...
                  timeVals <= params.baselineWindow(2);

    if ~any(baselineIdx)
        error('No samples found in baselineWindow.');
    end

    baselineVals = mean(dataIn(:,baselineIdx),2);
    dataIn = dataIn - baselineVals;
end

%% Time indices
subtractIdx = timeVals > params.subtractWindow(1) & ...
              timeVals <= params.subtractWindow(2);

scaleIdx = timeVals > params.scaleWindow(1) & ...
           timeVals <= params.scaleWindow(2);

if ~any(subtractIdx)
    error('No samples found in subtractWindow.');
end

if ~any(scaleIdx)
    error('No samples found in scaleWindow.');
end

%% ERP template
fullERP = mean(dataIn,1);

% Windowed template is what will actually be subtracted.
template = zeros(size(fullERP));
template(subtractIdx) = fullERP(subtractIdx);

% Use full ERP for scale estimation, not the zero-padded/windowed template.
templateFit = fullERP(scaleIdx);

denominator = sum(templateFit.^2);

if denominator == 0
    error('Template energy inside scaleWindow is zero. Cannot estimate scale.');
end

%% Trial-wise scaled template subtraction
numTrials = size(dataIn,1);

cleanedData = zeros(size(dataIn));
scaledTemplates = zeros(size(dataIn));
scaleFactors = zeros(numTrials,1);

for iTrial = 1:numTrials

    trial = dataIn(iTrial,:);
    trialFit = trial(scaleIdx);

    % Least-squares scale factor:
    % alpha = argmin ||trial - alpha*template||^2
    alpha = sum(trialFit .* templateFit) / denominator;

    % Bound scale factor to avoid unstable subtraction
    alpha = max(params.scaleBounds(1),min(params.scaleBounds(2),alpha));

    scaledTemplate = alpha * template;

    scaleFactors(iTrial) = alpha;
    scaledTemplates(iTrial,:) = scaledTemplate;
    cleanedData(iTrial,:) = trial - scaledTemplate;
end

%% Output
out = struct();
out.methodName = 'ERP scaled subtraction';
out.cleanedData = cleanedData;
out.template = template;
out.fullERP = fullERP;
out.scaledTemplates = scaledTemplates;
out.scaleFactors = scaleFactors;
out.rawMean = mean(dataIn,1);
out.cleanedMean = mean(cleanedData,1);
out.timeVals = timeVals;
out.params = params;

end