function out = ERPShiftedScaled(trialData,timeVals,params)
% ERPShiftedScaled
%
% ERP subtraction with trial-wise shift and scale correction.
%
% Important implementation choice:
%   - Use the full ERP template to estimate the best shift and scale.
%   - Subtract the shifted-scaled template only inside subtractWindow.
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

if ~isfield(params,'fitWindow') || isempty(params.fitWindow)
    params.fitWindow = [-0.01 0.03];
end

if ~isfield(params,'maxShiftMS') || isempty(params.maxShiftMS)
    params.maxShiftMS = 10;
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

fitIdx = timeVals > params.fitWindow(1) & ...
         timeVals <= params.fitWindow(2);

if ~any(subtractIdx)
    error('No samples found in subtractWindow.');
end

if ~any(fitIdx)
    error('No samples found in fitWindow.');
end

%% Shift candidates
dt = median(diff(timeVals));
maxShiftSamples = round((params.maxShiftMS/1000)/dt);
candidateShifts = -maxShiftSamples:maxShiftSamples;

%% ERP template
% Full ERP is used for shift and scale fitting.
fullERP = mean(dataIn,1);

% Windowed ERP is the nominal zero-shift subtraction template.
windowedERP = zeros(size(fullERP));
windowedERP(subtractIdx) = fullERP(subtractIdx);

numTrials = size(dataIn,1);
numTime = size(dataIn,2);

cleanedData = zeros(size(dataIn));
shiftedScaledTemplates = zeros(size(dataIn));
shiftSamples = zeros(numTrials,1);
shiftTimes = zeros(numTrials,1);
scaleFactors = zeros(numTrials,1);
fitError = zeros(numTrials,1);

%% Trial-wise shifted + scaled template subtraction
for iTrial = 1:numTrials

    trial = dataIn(iTrial,:);

    bestError = inf;
    bestShift = 0;
    bestScale = 1;
    bestTemplateForSubtraction = windowedERP;

    for iShift = 1:length(candidateShifts)

        shiftVal = candidateShifts(iShift);

        % Use full ERP for fitting.
        shiftedFullERP = shiftVectorNoWrap(fullERP,shiftVal);

        templateFit = shiftedFullERP(fitIdx);
        trialFit = trial(fitIdx);

        denom = sum(templateFit.^2);

        if denom == 0
            alpha = 0;
        else
            alpha = sum(trialFit .* templateFit) / denom;
        end

        % Bound scale factor
        alpha = max(params.scaleBounds(1),min(params.scaleBounds(2),alpha));

        candidateFit = alpha * templateFit;

        diffVals = trialFit - candidateFit;
        err = sum(diffVals.^2);

        if err < bestError
            bestError = err;
            bestShift = shiftVal;
            bestScale = alpha;

            % Subtract only inside the artifact/subtraction window.
            bestTemplateForSubtraction = zeros(1,numTime);
            bestTemplateForSubtraction(subtractIdx) = ...
                alpha * shiftedFullERP(subtractIdx);
        end
    end

    shiftSamples(iTrial) = bestShift;
    shiftTimes(iTrial) = bestShift * dt;
    scaleFactors(iTrial) = bestScale;
    fitError(iTrial) = bestError;

    shiftedScaledTemplates(iTrial,:) = bestTemplateForSubtraction;
    cleanedData(iTrial,:) = trial - bestTemplateForSubtraction;
end

%% Output
out = struct();
out.methodName = 'ERP shifted + scaled subtraction';
out.cleanedData = cleanedData;
out.template = windowedERP;
out.fullERP = fullERP;
out.shiftedScaledTemplates = shiftedScaledTemplates;
out.shiftSamples = shiftSamples;
out.shiftTimes = shiftTimes;
out.scaleFactors = scaleFactors;
out.fitError = fitError;
out.rawMean = mean(dataIn,1);
out.cleanedMean = mean(cleanedData,1);
out.timeVals = timeVals;
out.params = params;

end

function y = shiftVectorNoWrap(x,shiftSamples)
% Positive shift moves signal to the right.
% Negative shift moves signal to the left.
% No circular wrapping.

y = zeros(size(x));

if shiftSamples > 0
    y((shiftSamples+1):end) = x(1:(end-shiftSamples));

elseif shiftSamples < 0
    s = abs(shiftSamples);
    y(1:(end-s)) = x((s+1):end);

else
    y = x;
end

end