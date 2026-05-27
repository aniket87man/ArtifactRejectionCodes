function out = SMARTALite_nc(trialData, timeVals, opts)
% SMARTALite
%
% Simplified SMARTA-style artifact removal for trial-by-time data.
%
% Core ideas:
%   1) Build pulse-centered artifact segments.
%   2) For each segment, find K nearest segments across trials/pulses.
%   3) Use the median of nearest neighbors as the artifact template.
%   4) Subtract templates using overlap-add reconstruction.
%
% This is a SMARTA-inspired adaptation for ICMS trial data.
%
% INPUTS
%   trialData : [nTrials x nSamples]
%   timeVals  : [1 x nSamples] or [nSamples x 1], seconds
%   opts      : struct with fields
%
% Required opts:
%   opts.pulseTimes      : pulse times in seconds
%   opts.stiFreq         : stimulation frequency in Hz
%
% Optional opts:
%   opts.prePulse        : seconds before pulse, default 0.0005
%   opts.postPulseGap    : seconds before next pulse to stop segment, default 0.0005
%   opts.K               : number of nearest neighbors, default 20
%   opts.window          : FFT comparison window, default [0 0.4]
%   opts.computeFFT      : true/false, default true
%   opts.excludeSameTrial: true/false, default false
%
% OUTPUT
%   out.cleanedData      : cleaned trials, compatible with previous demos
%   out.cleanedTrials    : same as cleanedData
%   out.artifactEstimate : estimated artifact train
%   out.trialTemplates   : artifact templates per segment
%   out.diagnostics      : extra details

    %% --------------------------------------------------------------------
    % Defaults and checks
    % ---------------------------------------------------------------------

    if nargin < 3
        opts = struct();
    end

    if ~isfield(opts, 'pulseTimes') || isempty(opts.pulseTimes)
        error('opts.pulseTimes must be provided.');
    end

    if ~isfield(opts, 'stiFreq') || isempty(opts.stiFreq)
        error('opts.stiFreq must be provided.');
    end

    if ~isfield(opts, 'prePulse') || isempty(opts.prePulse)
        opts.prePulse = 0.0005; % 0.5 ms before pulse
    end

    if ~isfield(opts, 'postPulseGap') || isempty(opts.postPulseGap)
        opts.postPulseGap = 0.0005; % stop 0.5 ms before next pulse
    end

    if ~isfield(opts, 'K') || isempty(opts.K)
        opts.K = 20;
    end

    if ~isfield(opts, 'window') || isempty(opts.window)
        opts.window = [0 0.4];
    end

    if ~isfield(opts, 'computeFFT') || isempty(opts.computeFFT)
        opts.computeFFT = true;
    end

    if ~isfield(opts, 'excludeSameTrial') || isempty(opts.excludeSameTrial)
        opts.excludeSameTrial = false;
    end

    assert(isnumeric(trialData) && ismatrix(trialData), ...
        'trialData must be [nTrials x nSamples].');

    assert(isvector(timeVals), ...
        'timeVals must be a vector.');

    assert(size(trialData,2) == numel(timeVals), ...
        'trialData columns must match numel(timeVals).');

    timeVals = timeVals(:).';

    [nTrials, nSamples] = size(trialData);
    Fs = 1 / median(diff(timeVals));

    pulseTimes = opts.pulseTimes(:).';
    nPulsesOriginal = numel(pulseTimes);

    %% --------------------------------------------------------------------
    % Define segment around each pulse
    % ---------------------------------------------------------------------

    stPoint = -round(opts.prePulse * Fs);

    % Segment ends before the next expected pulse.
    % This avoids including the next pulse in the current segment.
    edPoint = round(((1 / opts.stiFreq) - opts.postPulseGap) * Fs);

    if edPoint <= 0
        error('Segment end is <= 0. Check opts.stiFreq and opts.postPulseGap.');
    end

    segLen = edPoint - stPoint + 1;

    if segLen <= 5
        error('Segment length is too short. Check prePulse, postPulseGap, and stiFreq.');
    end

    %% --------------------------------------------------------------------
    % Keep only pulse times whose full segment fits inside the trial
    % ---------------------------------------------------------------------

    validPulseMask = true(size(pulseTimes));

    for p = 1:nPulsesOriginal
        cIdx = nearestIndex(timeVals, pulseTimes(p));

        if (cIdx + stPoint < 1) || (cIdx + edPoint > nSamples)
            validPulseMask(p) = false;
        end
    end

    pulseTimes = pulseTimes(validPulseMask);
    nPulses = numel(pulseTimes);

    if nPulses == 0
        error('No valid pulseTimes remain after boundary checking.');
    end

    %% --------------------------------------------------------------------
    % Build all pulse segments across all trials and pulses
    % ---------------------------------------------------------------------

    nSeg = nTrials * nPulses;

    Xall = zeros(nSeg, segLen);
    segTrialIdx = zeros(nSeg,1);
    segPulseIdx = zeros(nSeg,1);
    segSampleIdx = zeros(nSeg, segLen);

    row = 0;

    for tr = 1:nTrials
        for p = 1:nPulses
            row = row + 1;

            cIdx = nearestIndex(timeVals, pulseTimes(p));
            idx = (cIdx + stPoint):(cIdx + edPoint);

            Xall(row,:) = trialData(tr, idx);
            segTrialIdx(row) = tr;
            segPulseIdx(row) = p;
            segSampleIdx(row,:) = idx;
        end
    end

    %% --------------------------------------------------------------------
    % Distance matrix
    % ---------------------------------------------------------------------

    % SMARTALite uses Euclidean distance on raw artifact segments.
    D = pairwiseEuclidean(Xall);

    %% --------------------------------------------------------------------
    % Estimate artifact template per segment
    % ---------------------------------------------------------------------

    Xhat = zeros(size(Xall));

    for i = 1:nSeg

        d = D(i,:);
        [~, order] = sort(d, 'ascend');

        % Exclude self
        order(order == i) = [];

        % Optionally exclude segments from same trial
        if opts.excludeSameTrial
            sameTrial = segTrialIdx(order) == segTrialIdx(i);
            order(sameTrial) = [];
        end

        if isempty(order)
            warning('No neighbor found for segment %d. Using all non-self segments.', i);
            order = setdiff(1:nSeg, i);
        end

        Kuse = min(opts.K, numel(order));
        neigh = order(1:Kuse);

        % Median template from nearest neighbors
        tpl = median(Xall(neigh,:), 1);

        % Taper segment edges for overlap-add reconstruction
        win = buildSegmentWindow( ...
            timeVals, pulseTimes, segPulseIdx(i), stPoint, edPoint, segLen);

        Xhat(i,:) = tpl .* win;
    end

    %% --------------------------------------------------------------------
    % Overlap-add reconstruction of artifact estimate and cleaned trials
    % ---------------------------------------------------------------------

    artifactEstimate = zeros(nTrials, nSamples);
    cleanedTrials = trialData;

    for i = 1:nSeg
        tr = segTrialIdx(i);
        idx = segSampleIdx(i,:);

        artifactEstimate(tr, idx) = artifactEstimate(tr, idx) + Xhat(i,:);
        cleanedTrials(tr, idx) = cleanedTrials(tr, idx) - Xhat(i,:);
    end

    %% --------------------------------------------------------------------
    % Diagnostics and FFT
    % ---------------------------------------------------------------------

    erp = mean(trialData, 1);
    cleanedERP = mean(cleanedTrials, 1);

    rmsRaw = sqrt(mean(trialData.^2, 1));
    rmsClean = sqrt(mean(cleanedTrials.^2, 1));

    stdRaw = std(trialData, 0, 1);
    stdClean = std(cleanedTrials, 0, 1);

    windowIdx = [];
    freqAxis = [];
    fftRawMean = [];
    fftCleanMean = [];

    if opts.computeFFT
        windowIdx = find(timeVals > opts.window(1) & timeVals < opts.window(2));

        if isempty(windowIdx)
            error('Selected FFT window has no samples.');
        end

        rawWin = trialData(:, windowIdx);
        cleanWin = cleanedTrials(:, windowIdx);

        nWin = numel(windowIdx);
        freqAxis = (0:nWin-1) * (Fs / nWin);

        fftRaw = abs(fft(rawWin, [], 2));
        fftClean = abs(fft(cleanWin, [], 2));

        fftRawMean = mean(fftRaw, 1);
        fftCleanMean = mean(fftClean, 1);
    end

    %% --------------------------------------------------------------------
    % Output
    % ---------------------------------------------------------------------

    out = struct();

    out.methodName = sprintf('SMARTALite_K%d', opts.K);

    out.rawTrials = trialData;
    out.cleanedTrials = cleanedTrials;

    % Compatibility with previous demos
    out.cleanedData = cleanedTrials;

    out.artifactEstimate = artifactEstimate;
    out.trialTemplates = Xhat;

    out.erp = erp;
    out.cleanedERP = cleanedERP;

    out.rawMean = erp;
    out.cleanedMean = cleanedERP;

    out.rmsRaw = rmsRaw;
    out.rmsClean = rmsClean;

    out.stdRaw = stdRaw;
    out.stdClean = stdClean;

    out.Fs = Fs;
    out.timeVals = timeVals;
    out.params = opts;

    out.window = opts.window;
    out.windowIdx = windowIdx;
    out.freqAxis = freqAxis;
    out.fftRawMean = fftRawMean;
    out.fftCleanMean = fftCleanMean;

    out.diagnostics = struct();
    out.diagnostics.originalPulseTimes = opts.pulseTimes;
    out.diagnostics.validPulseTimes = pulseTimes;
    out.diagnostics.nPulsesOriginal = nPulsesOriginal;
    out.diagnostics.nPulsesUsed = nPulses;
    out.diagnostics.segmentLength = segLen;
    out.diagnostics.stPoint = stPoint;
    out.diagnostics.edPoint = edPoint;
    out.diagnostics.segTrialIdx = segTrialIdx;
    out.diagnostics.segPulseIdx = segPulseIdx;
    out.diagnostics.segSampleIdx = segSampleIdx;
    out.diagnostics.segmentMatrix = Xall;
    out.diagnostics.artifactSegments = Xhat;
    out.diagnostics.artifactEstimate = artifactEstimate;
    out.diagnostics.K = opts.K;
end

%% ========================================================================
% Helper functions
% ========================================================================

function idx = nearestIndex(t, val)
    [~, idx] = min(abs(t - val));
end

function D = pairwiseEuclidean(X)
    G = sum(X.^2, 2);
    D2 = G + G' - 2*(X*X');
    D2(D2 < 0) = 0;
    D = sqrt(D2);
end

function win = buildSegmentWindow(timeVals, pulseTimes, pulseIdx, stPoint, edPoint, segLen)

    win = ones(1, segLen);

    nPulses = numel(pulseTimes);

    % Left overlap with previous pulse segment
    if pulseIdx ~= 1
        prevCenter = nearestIndex(timeVals, pulseTimes(pulseIdx-1));
        thisCenter = nearestIndex(timeVals, pulseTimes(pulseIdx));

        prevEnd = prevCenter + edPoint;
        thisStart = thisCenter + stPoint;

        overlap = prevEnd - thisStart + 1;

        if overlap > 1 && overlap < segLen
            taper = sin(linspace(0, pi/2, overlap)).^2;
            win(1:overlap) = taper;
        end
    end

    % Right overlap with next pulse segment
    if pulseIdx ~= nPulses
        thisCenter = nearestIndex(timeVals, pulseTimes(pulseIdx));
        nextCenter = nearestIndex(timeVals, pulseTimes(pulseIdx+1));

        thisEnd = thisCenter + edPoint;
        nextStart = nextCenter + stPoint;

        overlap = thisEnd - nextStart + 1;

        if overlap > 1 && overlap < segLen
            taper = cos(linspace(0, pi/2, overlap)).^2;
            win(end-overlap+1:end) = taper;
        end
    end
end