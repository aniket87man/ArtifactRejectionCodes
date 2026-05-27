function out = PulsewiseTemplate(trialData,timeVals,params)
% PulsewiseTemplate
%
% Pulse-wise stimulation artifact template subtraction.
%
% Main steps:
%   1. Detect pulse locations from the mean high-stimulation signal.
%   2. Estimate a local template around each detected pulse.
%   3. Combine local pulse templates into a full-length subtraction template.
%   4. Subtract the template only inside the selected subtraction window.
%
% Window convention:
%   [0 0.4] means timeVals > 0 and timeVals <= 0.4.
%   For 2000 Hz sampling, this selects 0.0005--0.4000 s.

if nargin < 3
    params = struct();
end

if ~isfield(params,'subtractWindow') || isempty(params.subtractWindow)
    params.subtractWindow = [timeVals(1)-eps timeVals(end)];
end

if ~isfield(params,'pulseSearchWindow') || isempty(params.pulseSearchWindow)
    params.pulseSearchWindow = [-0.01 0.32];
end

if ~isfield(params,'pulseWindowMS') || isempty(params.pulseWindowMS)
    params.pulseWindowMS = [-5 35];
end

if ~isfield(params,'minPulseDistanceMS') || isempty(params.minPulseDistanceMS)
    params.minPulseDistanceMS = 25;
end

if ~isfield(params,'expectedNumPulses')
    params.expectedNumPulses = [];
end

if ~isfield(params,'maxNumPulses') || isempty(params.maxNumPulses)
    params.maxNumPulses = 20;
end

if ~isfield(params,'thresholdMAD') || isempty(params.thresholdMAD)
    params.thresholdMAD = 6;
end

if ~isfield(params,'localBaselineMS') || isempty(params.localBaselineMS)
    params.localBaselineMS = 20;
end

if ~isfield(params,'templateStatistic') || isempty(params.templateStatistic)
    params.templateStatistic = 'median';
end

if ~isfield(params,'taperEdgeMS') || isempty(params.taperEdgeMS)
    params.taperEdgeMS = 2;
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

Fs = 1/median(diff(timeVals));

numTime = size(dataIn,2);

%% Time indices
subtractIdx = timeVals > params.subtractWindow(1) & ...
              timeVals <= params.subtractWindow(2);

searchIdx = find(timeVals > params.pulseSearchWindow(1) & ...
                 timeVals <= params.pulseSearchWindow(2));

if ~any(subtractIdx)
    error('No samples found in subtractWindow.');
end

if isempty(searchIdx)
    error('No samples found in pulseSearchWindow.');
end

%% Step 1: pulse detection from mean high-stimulation signal
rawMean = mean(dataIn,1);

localBaselineSamples = max(3,round((params.localBaselineMS/1000) * Fs));
localMean = movmean(rawMean,localBaselineSamples);

detectionScore = abs(rawMean - localMean);

scoreInWindow = detectionScore(searchIdx);
scoreMedian = median(scoreInWindow);
scoreMAD = median(abs(scoreInWindow - scoreMedian)) + eps;
threshold = scoreMedian + params.thresholdMAD * scoreMAD;

minPulseDistanceSamples = round((params.minPulseDistanceMS/1000) * Fs);

pulseSamples = detectLocalPeaksGreedy( ...
    detectionScore, ...
    searchIdx, ...
    threshold, ...
    minPulseDistanceSamples, ...
    params.expectedNumPulses, ...
    params.maxNumPulses);

if isempty(pulseSamples)
    error('No pulses detected. Try lowering params.thresholdMAD or setting params.expectedNumPulses.');
end

pulseSamples = sort(pulseSamples(:)');
pulseTimes = timeVals(pulseSamples);

%% Step 2: build local template around each pulse
pulseOffsets = round((params.pulseWindowMS(1)/1000) * Fs) : ...
               round((params.pulseWindowMS(2)/1000) * Fs);

numPulseSamples = length(pulseOffsets);
numPulses = length(pulseSamples);

pulseTemplates = zeros(numPulses,numPulseSamples);
fullTemplate = zeros(1,numTime);

edgeTaperSamples = round((params.taperEdgeMS/1000) * Fs);

for iPulse = 1:numPulses

    centerSample = pulseSamples(iPulse);
    templateSampleIdx = centerSample + pulseOffsets;

    validMask = templateSampleIdx >= 1 & templateSampleIdx <= numTime;
    validPositions = find(validMask);
    validIdx = templateSampleIdx(validMask);

    if isempty(validIdx)
        continue;
    end

    pulseSegments = dataIn(:,validIdx);

    switch lower(params.templateStatistic)
        case 'median'
            templateSegmentValid = median(pulseSegments,1);
        case 'mean'
            templateSegmentValid = mean(pulseSegments,1);
        otherwise
            error('params.templateStatistic must be ''median'' or ''mean''.');
    end

    templateSegment = zeros(1,numPulseSamples);
    templateSegment(validMask) = templateSegmentValid;

    templateSegment = applyEdgeTaper(templateSegment,edgeTaperSamples);

    pulseTemplates(iPulse,:) = templateSegment;

    % Add the local pulse template only inside the subtraction window.
    insideSubtractMask = subtractIdx(validIdx);

    addIdx = validIdx(insideSubtractMask);
    addPositions = validPositions(insideSubtractMask);

    fullTemplate(addIdx) = fullTemplate(addIdx) + templateSegment(addPositions);
end

%% Step 3: subtract pulse-wise template
cleanedData = dataIn - fullTemplate;

%% Output
out = struct();

out.methodName = 'Pulse-wise template subtraction';

out.cleanedData = cleanedData;
out.fullTemplate = fullTemplate;
out.pulseTemplates = pulseTemplates;

out.pulseSamples = pulseSamples;
out.pulseTimes = pulseTimes;
out.pulseOffsets = pulseOffsets;
out.pulseOffsetTimes = pulseOffsets / Fs;

out.detectionScore = detectionScore;
out.detectionThreshold = threshold;

out.rawMean = rawMean;
out.cleanedMean = mean(cleanedData,1);

out.timeVals = timeVals;
out.params = params;

end

function selectedPeaks = detectLocalPeaksGreedy(score,searchIdx,threshold,minDistanceSamples,expectedNumPulses,maxNumPulses)
% Detect local peaks above threshold, then choose the strongest peaks while
% enforcing a minimum separation.

scoreSearch = score(searchIdx);

localMaxMask = false(size(scoreSearch));

for i = 2:length(scoreSearch)-1
    if scoreSearch(i) >= scoreSearch(i-1) && ...
       scoreSearch(i) >  scoreSearch(i+1) && ...
       scoreSearch(i) >= threshold
        localMaxMask(i) = true;
    end
end

candidateSamples = searchIdx(localMaxMask);

% If too few candidates are found and expectedNumPulses is specified,
% relax the threshold and use all local maxima.
if ~isempty(expectedNumPulses) && length(candidateSamples) < expectedNumPulses

    localMaxMask = false(size(scoreSearch));

    for i = 2:length(scoreSearch)-1
        if scoreSearch(i) >= scoreSearch(i-1) && scoreSearch(i) > scoreSearch(i+1)
            localMaxMask(i) = true;
        end
    end

    candidateSamples = searchIdx(localMaxMask);
end

if isempty(candidateSamples)
    selectedPeaks = [];
    return;
end

candidateScores = score(candidateSamples);
[~,order] = sort(candidateScores,'descend');

selectedPeaks = [];

for i = 1:length(order)

    candidate = candidateSamples(order(i));

    if isempty(selectedPeaks)
        selectedPeaks = candidate;
    else
        if all(abs(candidate - selectedPeaks) >= minDistanceSamples)
            selectedPeaks(end+1) = candidate; %#ok<AGROW>
        end
    end

    if ~isempty(expectedNumPulses)
        if length(selectedPeaks) >= expectedNumPulses
            break;
        end
    else
        if length(selectedPeaks) >= maxNumPulses
            break;
        end
    end
end

selectedPeaks = sort(selectedPeaks);

end

function y = applyEdgeTaper(x,edgeSamples)
% Apply a raised-cosine taper to both edges of a local template.

y = x;

if edgeSamples <= 0
    return;
end

n = length(x);

edgeSamples = min(edgeSamples,floor(n/2));

if edgeSamples <= 0
    return;
end

t = linspace(0,pi/2,edgeSamples);
leftTaper = sin(t).^2;
rightTaper = fliplr(leftTaper);

y(1:edgeSamples) = y(1:edgeSamples) .* leftTaper;
y(end-edgeSamples+1:end) = y(end-edgeSamples+1:end) .* rightTaper;

end