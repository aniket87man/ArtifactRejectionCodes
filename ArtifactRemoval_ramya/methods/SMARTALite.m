function out = SMARTALite(trialData, timeVals, opts)
% SMARTALite
%
% Simplified SMARTA-style artifact removal for trial-by-time data.
%
% Core ideas:
%   1) build inter-pulse segments
%   2) for each segment, find K nearest segments across all trials
%   3) use median of neighbors as template
%   4) reconstruct with overlap-add tapering
%
% INPUTS
%   trialData : [nTrials x nSamples]
%   timeVals  : [1 x nSamples] or [nSamples x 1]
%   opts      : struct with fields
%       .pulseTimes      = pulse times in seconds
%       .stiFreq         = stimulation frequency in Hz
%       .prePulse        = seconds before pulse, default 0.0005
%       .K               = number of neighbors, default 20
%       .window          = FFT comparison window, default [0 0.4]
%       .computeFFT      = true/false, default true
%
% OUTPUT
%   out : struct with standard fields

    if nargin < 3
        opts = struct();
    end
    if ~isfield(opts, 'pulseTimes')
        error('opts.pulseTimes must be provided.');
    end
    if ~isfield(opts, 'stiFreq')
        error('opts.stiFreq must be provided.');
    end
    if ~isfield(opts, 'prePulse')
        opts.prePulse = 0.0005; % 0.5 ms
    end
    if ~isfield(opts, 'K')
        opts.K = 5;
    end
    if ~isfield(opts, 'window')
        opts.window = [0 0.4];
    end
    if ~isfield(opts, 'computeFFT')
        opts.computeFFT = true;
    end
    if ~isfield(opts, 'excludeSameTrial')
        opts.excludeSameTrial = true;
    end

    assert(isnumeric(trialData) && ismatrix(trialData), ...
        'trialData must be [nTrials x nSamples].');
    assert(isvector(timeVals), 'timeVals must be a vector.');
    assert(size(trialData,2) == numel(timeVals), ...
        'trialData columns must match numel(timeVals).');

    timeVals = timeVals(:).';
    [nTrials, nSamples] = size(trialData);
    Fs = 1 / median(diff(timeVals));

    pulseTimes = opts.pulseTimes(:).';
    nPulses = numel(pulseTimes);

    stPoint = -round(opts.prePulse * Fs);
    edPoint = round((1 / opts.stiFreq) * Fs) - 1;
    segLen = edPoint - stPoint + 1;
    midPoint = round(segLen / 2);

    % Keep only pulse times whose full segment fits in trial
    validPulseMask = true(size(pulseTimes));
    for p = 1:nPulses
        cIdx = nearestIndex(timeVals, pulseTimes(p));
        if (cIdx + stPoint < 1) || (cIdx + edPoint > nSamples)
            validPulseMask(p) = false;
        end
    end
    pulseTimes = pulseTimes(validPulseMask);
    nPulses = numel(pulseTimes);

    % Build all segments across all trials and all pulses
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

    % Pairwise Euclidean distances between segments
    % pdist2 may not be available in all MATLAB installs, so use manual formula
    D = pairwiseEuclidean(Xall);

    % Artifact estimate per segment
    Xhat = zeros(size(Xall));

    for i = 1:nSeg
        d = D(i,:);
        [~, order] = sort(d, 'ascend');

        % Exclude self
        order(order == i) = [];

        % Optionally exclude segments from the same trial
        if opts.excludeSameTrial
            sameTrialMask = segTrialIdx(order) == segTrialIdx(i);
            order(sameTrialMask) = [];
        end


        % Safety fallback
        if isempty(order)
            warning('No eligible neighbors found after excluding same-trial segments. Using all non-self segments.');
            order = 1:nSeg;
            order(order == i) = [];
        end
        
        Kuse = min(opts.K, numel(order));
        neigh = order(1:Kuse);

        % median template from nearest neighbors
        tpl = median(Xall(neigh,:), 1);

        % overlap-add tapering similar to SMARTA code
        p = segPulseIdx(i);
        tr = segTrialIdx(i);

        win = ones(1, segLen);

        % Left overlap
        if p ~= 1
            prevCenter = nearestIndex(timeVals, pulseTimes(p-1));
            thisCenter = nearestIndex(timeVals, pulseTimes(p));
            gap = edPoint - (thisCenter - prevCenter + stPoint);
            if gap > 0 && gap < segLen
                win(1:gap) = sin(pi*(0:gap-1)/(2*gap)).^2;
            end
        end

        % Right overlap
        if p ~= nPulses
            thisCenter = nearestIndex(timeVals, pulseTimes(p));
            nextCenter = nearestIndex(timeVals, pulseTimes(p+1));
            gap = edPoint - (nextCenter - thisCenter + stPoint);
            if gap > 0 && gap < segLen
                win(end-gap+1:end) = cos(pi*(1:gap)/(2*gap)).^2;
            end
        end

        Xhat(i,:) = tpl .* win;
    end

    % Overlap-add reconstruction of artifact estimate
    sa = zeros(nTrials, nSamples);
    z = trialData;

    for i = 1:nSeg
        tr = segTrialIdx(i);
        idx = segSampleIdx(i,:);

        sa(tr, idx) = sa(tr, idx) + Xhat(i,:);
        z(tr, idx) = z(tr, idx) - Xhat(i,:);
    end

    cleanedTrials = z;

    % Diagnostics
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
        windowIdx = find(timeVals > opts.window(1) & timeVals <= opts.window(2));
        assert(~isempty(windowIdx), 'Selected FFT window has no samples.');

        rawWin = trialData(:, windowIdx);
        cleanWin = cleanedTrials(:, windowIdx);

        nWin = numel(windowIdx);
        freqAxis = (0:nWin-1) * (Fs / nWin);

        fftRawMean = log10(mean(abs(fft(rawWin'))'));
        fftCleanMean = log10(mean(abs(fft(cleanWin'))'));
    end

    out = struct();
    out.methodName = sprintf('SMARTALite_K%d', opts.K);
    out.rawTrials = trialData;
    out.cleanedTrials = cleanedTrials;
    out.cleanedData = cleanedTrials;
    out.artifactEstimate = sa;
    out.trialTemplates = Xhat;
    out.params = opts;
    out.timeVals = timeVals;
    out.erp = erp;
    out.cleanedERP = cleanedERP;
    out.rmsRaw = rmsRaw;
    out.rmsClean = rmsClean;
    out.stdRaw = stdRaw;
    out.stdClean = stdClean;
    out.Fs = Fs;
    out.window = opts.window;
    out.windowIdx = windowIdx;
    out.freqAxis = freqAxis;
    out.fftRawMean = fftRawMean;
    out.fftCleanMean = fftCleanMean;

    out.diagnostics = struct();
    out.diagnostics.pulseTimes = pulseTimes;
    out.diagnostics.segmentLength = segLen;
    out.diagnostics.segTrialIdx = segTrialIdx;
    out.diagnostics.segPulseIdx = segPulseIdx;
    out.diagnostics.segSampleIdx = segSampleIdx;
    out.diagnostics.segmentMatrix = Xall;
    out.diagnostics.artifactSegments = Xhat;
    out.diagnostics.artifactEstimate = sa;
    out.diagnostics.K = opts.K;
    out.diagnostics.midPoint = midPoint;
end

function idx = nearestIndex(t, val)
    [~, idx] = min(abs(t - val));
end

function D = pairwiseEuclidean(X)
    G = sum(X.^2, 2);
    D2 = G + G' - 2*(X*X');
    D2(D2 < 0) = 0;
    D = sqrt(D2);
end