function out = ERPShifted(dataIn,timeVals,params)
% ERPShifted
% Trial-wise shifted ERP/template subtraction.
%
% Important implementation choice:
%   - Use the full ERP template to estimate the best temporal shift.
%   - Subtract the shifted template only inside params.subtractWindow.
%
% This avoids modifying the full trial, but still estimates the shift from
% the real ERP waveform instead of a window-masked template.

%% Defaults
if nargin < 3
    params = struct();
end

if ~isfield(params,'subtractWindow') || isempty(params.subtractWindow)
    params.subtractWindow = [timeVals(1) timeVals(end)];
end

if ~isfield(params,'alignWindow') || isempty(params.alignWindow)
    params.alignWindow = params.subtractWindow;
end

if ~isfield(params,'maxShiftMS') || isempty(params.maxShiftMS)
    params.maxShiftMS = 10;
end

if ~isfield(params,'doBaselineCorrection') || isempty(params.doBaselineCorrection)
    params.doBaselineCorrection = false;
end

%% Basic checks
if size(dataIn,2) ~= length(timeVals)
    error('Number of columns in dataIn must match length(timeVals).');
end

dataWork = dataIn;

%% Optional baseline correction
if params.doBaselineCorrection
    baselineIdx = timeVals < 0;
    baselineVals = mean(dataWork(:,baselineIdx),2);
    dataWork = dataWork - baselineVals;
end

%% Time indices
subtractIdx = timeVals > params.subtractWindow(1) & ...
              timeVals <= params.subtractWindow(2);

alignIdx = timeVals > params.alignWindow(1) & ...
           timeVals <= params.alignWindow(2);

if ~any(subtractIdx)
    error('subtractWindow does not contain any samples.');
end

if ~any(alignIdx)
    error('alignWindow does not contain any samples.');
end

%% ERP template
% Full ERP is used for shift estimation.
fullERP = mean(dataWork,1);

% Windowed ERP is the nominal component subtracted in the zero-shift case.
windowedERP = zeros(size(fullERP));
windowedERP(subtractIdx) = fullERP(subtractIdx);

%% Candidate shifts
dt = median(diff(timeVals));
maxShiftSamples = round((params.maxShiftMS/1000)/dt);
candidateShifts = -maxShiftSamples:maxShiftSamples;

nTrials = size(dataWork,1);
nTime = size(dataWork,2);

cleanedData = zeros(size(dataWork));
shiftedTemplates = zeros(size(dataWork));
shiftSamples = zeros(nTrials,1);
shiftTimes = zeros(nTrials,1);
alignmentError = zeros(nTrials,1);

%% Trial-wise shifted template subtraction
for iTrial = 1:nTrials

    trial = dataWork(iTrial,:);

    bestError = inf;
    bestShift = 0;
    bestTemplateForSubtraction = windowedERP;

    for iShift = 1:length(candidateShifts)

        shiftVal = candidateShifts(iShift);

        % Corrected part:
        % Use the full ERP for alignment/shift estimation.
        shiftedFullERP = shiftVectorNoWrap(fullERP,shiftVal);

        diffVals = trial(alignIdx) - shiftedFullERP(alignIdx);
        err = sum(diffVals.^2);

        if err < bestError
            bestError = err;
            bestShift = shiftVal;

            % But subtract only inside the subtraction window.
            bestTemplateForSubtraction = zeros(1,nTime);
            bestTemplateForSubtraction(subtractIdx) = shiftedFullERP(subtractIdx);
        end
    end

    shiftSamples(iTrial) = bestShift;
    shiftTimes(iTrial) = bestShift * dt;
    alignmentError(iTrial) = bestError;

    shiftedTemplates(iTrial,:) = bestTemplateForSubtraction;
    cleanedData(iTrial,:) = trial - bestTemplateForSubtraction;
end

%% Output
out = struct();
out.cleanedData = cleanedData;
out.template = fullERP;
out.windowedTemplate = windowedERP;
out.shiftedTemplates = shiftedTemplates;
out.shiftSamples = shiftSamples;
out.shiftTimes = shiftTimes;
out.alignmentError = alignmentError;
out.params = params;
out.methodName = 'ERP shifted subtraction';
out.fullERP = fullERP;
out.rawMean = mean(dataWork,1);
out.cleanedMean = mean(cleanedData,1);
out.timeVals = timeVals;

end

%% Helper function
function y = shiftVectorNoWrap(x,shiftSamples)
% Shift vector without circular wrapping.
% Positive shift moves the signal to the right.
% Negative shift moves the signal to the left.

n = length(x);
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