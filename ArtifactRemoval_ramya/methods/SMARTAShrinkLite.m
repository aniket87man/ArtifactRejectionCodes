function out = SMARTAShrinkLite(trialData,timeVals,opts)
% SMARTAShrinkLite
%
% SMARTA-inspired stimulation-cycle artifact removal with SVD-shrinkage
% similarity.
%
% This is a modified version of SMARTALite:
%   1. Segment each trial into full stimulation-cycle segments.
%   2. Build a segment matrix across all trials and pulses.
%   3. Use SVD shrinkage only to denoise the segment matrix for distance/KNN.
%   4. Estimate each artifact template from original nearest-neighbor segments.
%   5. Subtract the template from the original trial data.
%
% Important:
%   K = number of nearest-neighbor segments for template estimation.
%   K is NOT number of PCA components removed.
%
% Required opts:
%   opts.pulseTimes : pulse/cycle start times in seconds
%   opts.stiFreq    : stimulation frequency in Hz
%
% Recommended for current ICMS data:
%   opts.pulseTimes = 0 + (0:7)*0.05;
%   opts.stiFreq = 20;
%   opts.prePulse = 0.0005;
%   opts.K = 3;
%   opts.maxRank = 10;
%   opts.shrinkStrength = 1.0;
%   opts.hpCutoff = 0 or 150;

%% ========================================================================
% Defaults and checks
% ========================================================================

if nargin < 3
    opts = struct();
end

if ~isfield(opts,'pulseTimes') || isempty(opts.pulseTimes)
    error('opts.pulseTimes must be provided.');
end

if ~isfield(opts,'stiFreq') || isempty(opts.stiFreq)
    error('opts.stiFreq must be provided.');
end

if ~isfield(opts,'prePulse') || isempty(opts.prePulse)
    opts.prePulse = 0.0005;
end

if ~isfield(opts,'K') || isempty(opts.K)
    opts.K = 3;
end

if ~isfield(opts,'hpCutoff') || isempty(opts.hpCutoff)
    opts.hpCutoff = 0; % 0 means no high-pass before shrinkage
end

if ~isfield(opts,'maxRank') || isempty(opts.maxRank)
    opts.maxRank = 10;
end

if ~isfield(opts,'minRank') || isempty(opts.minRank)
    opts.minRank = 2;
end

if ~isfield(opts,'shrinkStrength') || isempty(opts.shrinkStrength)
    opts.shrinkStrength = 1.0;
end

if ~isfield(opts,'useSoftShrinkage') || isempty(opts.useSoftShrinkage)
    opts.useSoftShrinkage = false;
end

if ~isfield(opts,'excludeSelf') || isempty(opts.excludeSelf)
    opts.excludeSelf = true;
end

if ~isfield(opts,'excludeSameTrial') || isempty(opts.excludeSameTrial)
    opts.excludeSameTrial = false;
end

if ~isfield(opts,'computeFFT') || isempty(opts.computeFFT)
    opts.computeFFT = false;
end

if ~isfield(opts,'window') || isempty(opts.window)
    opts.window = [0 0.4];
end

if size(trialData,2) ~= length(timeVals)
    error('trialData must be trials x time, and length(timeVals) must match size(trialData,2).');
end

timeVals = timeVals(:).';

[nTrials,nSamples] = size(trialData);
Fs = 1/median(diff(timeVals));

pulseTimes = opts.pulseTimes(:).';
nPulsesOriginal = length(pulseTimes);

%% ========================================================================
% Segment definition
% ========================================================================

stPoint = -round(opts.prePulse * Fs);

% Old successful SMARTALite behavior:
% full stimulation-cycle segment from prePulse before pulse to 1/stiFreq after.
edPoint = round((1/opts.stiFreq) * Fs);

segLen = edPoint - stPoint + 1;

if segLen <= 5
    error('Segment length is too short. Check opts.prePulse and opts.stiFreq.');
end

%% ========================================================================
% Keep only pulseTimes whose full segment fits in the trial
% ========================================================================

validPulseMask = true(size(pulseTimes));

for iPulse = 1:nPulsesOriginal
    centerIdx = nearestIndex(timeVals,pulseTimes(iPulse));

    if centerIdx + stPoint < 1 || centerIdx + edPoint > nSamples
        validPulseMask(iPulse) = false;
    end
end

pulseTimes = pulseTimes(validPulseMask);
nPulses = length(pulseTimes);

if nPulses == 0
    error('No pulseTimes remain after boundary checking.');
end

%% ========================================================================
% Build segment matrix
% ========================================================================

nSegments = nTrials * nPulses;

X = zeros(nSegments,segLen);
segmentTrial = zeros(nSegments,1);
segmentPulse = zeros(nSegments,1);
segmentSamples = zeros(nSegments,segLen);

rowCounter = 0;

for iTrial = 1:nTrials
    for iPulse = 1:nPulses

        rowCounter = rowCounter + 1;

        centerIdx = nearestIndex(timeVals,pulseTimes(iPulse));
        idx = (centerIdx + stPoint):(centerIdx + edPoint);

        X(rowCounter,:) = trialData(iTrial,idx);
        segmentTrial(rowCounter) = iTrial;
        segmentPulse(rowCounter) = iPulse;
        segmentSamples(rowCounter,:) = idx;
    end
end

%% ========================================================================
% Build denoised feature matrix for KNN
% ========================================================================

Xfeature = X;

% Optional high-pass filtering for similarity only
if opts.hpCutoff > 0
    if opts.hpCutoff >= Fs/2
        warning('hpCutoff >= Nyquist. Skipping high-pass filtering.');
    else
        try
            [b,a] = butter(3,opts.hpCutoff/(Fs/2),'high');
            Xfeature = filtfilt(b,a,Xfeature')';
        catch ME
            warning('High-pass filtering failed: %s. Continuing without high-pass.',ME.message);
            Xfeature = X;
        end
    end
end

% Remove common mean waveform before SVD.
% This does not remove anything from final data; it only improves neighbor search.
segmentMean = mean(Xfeature,1);
Xcentered = Xfeature - segmentMean;

% SVD shrinkage
[U,S,V] = svd(Xcentered,'econ');
singVals = diag(S);

if isempty(singVals)
    error('SVD failed: no singular values.');
end

% Gavish-Donoho-style simple threshold estimate.
% This is not exact eOptShrink. It is a stable modified shrinkage heuristic.
m = size(Xcentered,1);
n = size(Xcentered,2);
beta = min(m,n) / max(m,n);

omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
tau = opts.shrinkStrength * omega * median(singVals);

rankEstimated = sum(singVals > tau);

rankUse = max(opts.minRank,rankEstimated);
rankUse = min(rankUse,opts.maxRank);
rankUse = min(rankUse,length(singVals));

if rankUse < 1
    rankUse = min(opts.maxRank,length(singVals));
end

if opts.useSoftShrinkage
    sKeep = max(singVals(1:rankUse) - tau,0);

    % If soft shrinkage kills all values, fall back to hard rank.
    if all(sKeep == 0)
        sKeep = singVals(1:rankUse);
    end
else
    sKeep = singVals(1:rankUse);
end

Xdenoised = U(:,1:rankUse) * diag(sKeep) * V(:,1:rankUse)';

% Add mean back. Pairwise distances would not require this, but it keeps
% feature scale interpretable.
Xdenoised = Xdenoised + segmentMean;

%% ========================================================================
% KNN using denoised features, template from original segments
% ========================================================================

D = pairwiseEuclidean(Xdenoised);

Xhat = zeros(size(X));

for iSegment = 1:nSegments

    d = D(iSegment,:);
    [~,order] = sort(d,'ascend');

    if opts.excludeSelf
        order(order == iSegment) = [];
    end

    if opts.excludeSameTrial
        sameTrial = segmentTrial(order) == segmentTrial(iSegment);
        order(sameTrial) = [];
    end

    if isempty(order)
        order = setdiff(1:nSegments,iSegment);
    end

    Kuse = min(opts.K,length(order));
    neighborIdx = order(1:Kuse);

    templateSegment = median(X(neighborIdx,:),1);

    win = buildSegmentWindow(timeVals,pulseTimes,segmentPulse(iSegment),stPoint,edPoint,segLen);

    Xhat(iSegment,:) = templateSegment .* win;
end

%% ========================================================================
% Reconstruct artifact train and cleaned trials
% ========================================================================

sa = zeros(nTrials,nSamples);
cleanedTrials = trialData;

for iSegment = 1:nSegments
    iTrial = segmentTrial(iSegment);
    idx = segmentSamples(iSegment,:);

    sa(iTrial,idx) = sa(iTrial,idx) + Xhat(iSegment,:);
    cleanedTrials(iTrial,idx) = cleanedTrials(iTrial,idx) - Xhat(iSegment,:);
end

%% ========================================================================
% Optional FFT summaries
% ========================================================================

freqAxis = [];
fftRawMean = [];
fftCleanMean = [];
windowIdx = [];

if opts.computeFFT
    windowIdx = find(timeVals > opts.window(1) & timeVals < opts.window(2));

    if ~isempty(windowIdx)
        nWin = length(windowIdx);
        freqAxis = (0:nWin-1) * (Fs/nWin);

        fftRaw = abs(fft(trialData(:,windowIdx),[],2));
        fftClean = abs(fft(cleanedTrials(:,windowIdx),[],2));

        fftRawMean = mean(fftRaw,1);
        fftCleanMean = mean(fftClean,1);
    end
end

%% ========================================================================
% Output
% ========================================================================

out = struct();

out.methodName = sprintf('SMARTAShrinkLite_K%d_rank%d',opts.K,rankUse);

out.rawTrials = trialData;
out.cleanedTrials = cleanedTrials;

% Compatibility with demo pipeline
out.cleanedData = cleanedTrials;

out.artifactEstimate = sa;
out.trialTemplates = Xhat;

out.rawMean = mean(trialData,1);
out.cleanedMean = mean(cleanedTrials,1);

out.erp = out.rawMean;
out.cleanedERP = out.cleanedMean;

out.timeVals = timeVals;
out.Fs = Fs;
out.params = opts;

out.freqAxis = freqAxis;
out.fftRawMean = fftRawMean;
out.fftCleanMean = fftCleanMean;
out.windowIdx = windowIdx;

out.diagnostics = struct();
out.diagnostics.originalPulseTimes = opts.pulseTimes;
out.diagnostics.validPulseTimes = pulseTimes;
out.diagnostics.nPulsesOriginal = nPulsesOriginal;
out.diagnostics.nPulsesUsed = nPulses;
out.diagnostics.stPoint = stPoint;
out.diagnostics.edPoint = edPoint;
out.diagnostics.segLen = segLen;
out.diagnostics.segmentMatrix = X;
out.diagnostics.featureMatrix = Xfeature;
out.diagnostics.denoisedFeatureMatrix = Xdenoised;
out.diagnostics.singularValues = singVals;
out.diagnostics.threshold = tau;
out.diagnostics.rankEstimated = rankEstimated;
out.diagnostics.rankUse = rankUse;
out.diagnostics.segmentTrial = segmentTrial;
out.diagnostics.segmentPulse = segmentPulse;
out.diagnostics.segmentSamples = segmentSamples;
out.diagnostics.K = opts.K;

end

%% ========================================================================
% Helper functions
% ========================================================================

function idx = nearestIndex(t,val)
    [~,idx] = min(abs(t-val));
end

function D = pairwiseEuclidean(X)
    G = sum(X.^2,2);
    D2 = G + G' - 2*(X*X');
    D2(D2 < 0) = 0;
    D = sqrt(D2);
end

function win = buildSegmentWindow(timeVals,pulseTimes,pulseIdx,stPoint,edPoint,segLen)

    win = ones(1,segLen);
    nPulses = length(pulseTimes);

    % Left overlap with previous segment
    if pulseIdx > 1
        prevCenter = nearestIndex(timeVals,pulseTimes(pulseIdx-1));
        thisCenter = nearestIndex(timeVals,pulseTimes(pulseIdx));

        prevEnd = prevCenter + edPoint;
        thisStart = thisCenter + stPoint;

        overlap = prevEnd - thisStart + 1;

        if overlap > 1 && overlap < segLen
            taper = sin(linspace(0,pi/2,overlap)).^2;
            win(1:overlap) = taper;
        end
    end

    % Right overlap with next segment
    if pulseIdx < nPulses
        thisCenter = nearestIndex(timeVals,pulseTimes(pulseIdx));
        nextCenter = nearestIndex(timeVals,pulseTimes(pulseIdx+1));

        thisEnd = thisCenter + edPoint;
        nextStart = nextCenter + stPoint;

        overlap = thisEnd - nextStart + 1;

        if overlap > 1 && overlap < segLen
            taper = cos(linspace(0,pi/2,overlap)).^2;
            win(end-overlap+1:end) = taper;
        end
    end
end