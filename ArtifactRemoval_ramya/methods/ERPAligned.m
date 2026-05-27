function out = ERPAligned(trialData,timeVals,params)
% ERPAligned
%
% ERP template subtraction after trial-wise alignment.
%
% Main idea:
% 1. Detect the first stimulation pulse timing in each trial.
% 2. Align all trials to a common pulse timing.
% 3. Compute ERP/template from aligned trials.
% 4. Shift the full aligned template back to each trial's original timing.
% 5. Subtract only inside params.subtractWindow.
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

if ~isfield(params,'alignWindow') || isempty(params.alignWindow)
    params.alignWindow = [-0.01 0.03];
end

if ~isfield(params,'maxShiftMS') || isempty(params.maxShiftMS)
    params.maxShiftMS = 10;
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
Fs = 1/median(diff(timeVals));
maxShiftSamples = round((params.maxShiftMS/1000) * Fs);

subtractIdx = timeVals > params.subtractWindow(1) & ...
              timeVals <= params.subtractWindow(2);

alignIdx = find(timeVals > params.alignWindow(1) & ...
                timeVals <= params.alignWindow(2));

if ~any(subtractIdx)
    error('No samples found in subtractWindow.');
end

if isempty(alignIdx)
    error('No samples found in alignWindow.');
end

numTrials = size(dataIn,1);
numTime = size(dataIn,2);

%% Step 1: detect first-pulse location in each trial
peakSamples = zeros(numTrials,1);

for iTrial = 1:numTrials

    trialSegment = dataIn(iTrial,alignIdx);

    % Remove local median so peak detection focuses on sharp deflection.
    trialSegment = trialSegment - median(trialSegment);

    % Detect largest absolute deflection inside alignWindow.
    [~,localPeakIdx] = max(abs(trialSegment));

    peakSamples(iTrial) = alignIdx(localPeakIdx);
end

% Common alignment target: median peak location across trials.
targetPeakSample = round(median(peakSamples));

%% Step 2: align all trials
alignedTrials = zeros(size(dataIn));
shiftSamples = zeros(numTrials,1);
wasClipped = false(numTrials,1);

for iTrial = 1:numTrials

    rawShift = targetPeakSample - peakSamples(iTrial);

    % Limit extreme shifts.
    clippedShift = max(-maxShiftSamples,min(maxShiftSamples,rawShift));

    if clippedShift ~= rawShift
        wasClipped(iTrial) = true;
    end

    shiftSamples(iTrial) = clippedShift;
    alignedTrials(iTrial,:) = shiftVectorNoWrap(dataIn(iTrial,:),clippedShift);
end

%% Step 3: compute full ERP/template from aligned trials
fullAlignedERP = mean(alignedTrials,1);

% Windowed aligned template is useful for diagnostics only.
template = zeros(1,numTime);
template(subtractIdx) = fullAlignedERP(subtractIdx);

%% Step 4: shift full template back and subtract only inside window
trialTemplates = zeros(size(dataIn));
cleanedData = zeros(size(dataIn));

for iTrial = 1:numTrials

    % Since trial was aligned using +shiftSamples(iTrial),
    % move template back using -shiftSamples(iTrial).
    shiftedBackFullERP = shiftVectorNoWrap(fullAlignedERP,-shiftSamples(iTrial));

    % Important correction:
    % Subtract only inside subtractWindow after shifting back.
    trialTemplate = zeros(1,numTime);
    trialTemplate(subtractIdx) = shiftedBackFullERP(subtractIdx);

    trialTemplates(iTrial,:) = trialTemplate;
    cleanedData(iTrial,:) = dataIn(iTrial,:) - trialTemplate;
end

%% Output
out = struct();

out.methodName = 'ERP aligned subtraction';

out.cleanedData = cleanedData;
out.template = template;
out.fullAlignedERP = fullAlignedERP;

out.alignedTrials = alignedTrials;
out.trialTemplates = trialTemplates;

out.shiftSamples = shiftSamples;
out.shiftTimes = shiftSamples / Fs;
out.peakSamples = peakSamples;
out.peakTimes = timeVals(peakSamples);
out.targetPeakSample = targetPeakSample;
out.targetPeakTime = timeVals(targetPeakSample);
out.wasClipped = wasClipped;

out.rawMean = mean(dataIn,1);
out.cleanedMean = mean(cleanedData,1);

out.timeVals = timeVals;
out.params = params;

end

function y = shiftVectorNoWrap(x,shiftSamples)
% shiftVectorNoWrap
%
% Positive shift moves signal to the right.
% Negative shift moves signal to the left.
% No circular wrapping is used. Empty samples are filled with zero.

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