function out = SMARTAFull(trialData, timeVals, opts)
% SMARTAFull
%
% Fuller SMARTA-style artifact removal for trial-by-time data.
%
% Differences from SMARTALite:
%   - uses high-passed inter-pulse segments for similarity
%   - uses optimal_shrinkage_color_fast for denoising before distance calc
%
% Keeps:
%   - manual pulse times
%   - fixed K nearest neighbors
%   - median template from ORIGINAL segments
%   - overlap-add tapering
%
% INPUTS
%   trialData : [nTrials x nSamples]
%   timeVals  : [1 x nSamples] or [nSamples x 1]
%   opts      : struct with fields
%       .pulseTimes     = pulse times in seconds
%       .stiFreq        = stimulation frequency in Hz
%       .prePulse       = seconds before pulse, default 0.0005
%       .K              = number of neighbors, default 10
%       .window         = FFT comparison window, default [0 0.4]
%       .computeFFT     = true/false, default true
%       .hpCutoff       = high-pass cutoff for similarity, default 300
%       .shrinkLoss     = default 'fro'
%       .shrinkKL       = default 10
%       .shrinkKH       = default 15
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
        opts.prePulse = 0.0005;
    end
    if ~isfield(opts, 'K')
        opts.K = 10;
    end
    if ~isfield(opts, 'window')
        opts.window = [0 0.4];
    end
    if ~isfield(opts, 'computeFFT')
        opts.computeFFT = true;
    end
    if ~isfield(opts, 'hpCutoff')
        opts.hpCutoff = 300;
    end
    if ~isfield(opts, 'shrinkLoss')
        opts.shrinkLoss = 'fro';
    end
    if ~isfield(opts, 'shrinkKL')
        opts.shrinkKL = 10;
    end
    if ~isfield(opts, 'shrinkKH')
        opts.shrinkKH = 15;
    end

    assert(isnumeric(trialData) && ndims(trialData) == 2, ...
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
    edPoint = round((1 / opts.stiFreq) * Fs);
    segLen = edPoint - stPoint + 1;

    % Keep only pulse times whose full segment fits
    validPulseMask = true(size(pulseTimes));
    for p = 1:nPulses
        cIdx = nearestIndex(timeVals, pulseTimes(p));
        if (cIdx + stPoint < 1) || (cIdx + edPoint > nSamples)
            validPulseMask(p) = false;
        end
    end
    pulseTimes = pulseTimes(validPulseMask);
    nPulses = numel(pulseTimes);

    % Build all segments
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

    % High-pass segments for similarity
    [bHp, aHp] = butter(3, opts.hpCutoff / (Fs/2), 'high');
    Xhp = filtfilt(bHp, aHp, Xall')';
    
    % Optimal shrinkage denoising for similarity
    usedShrinkage = true;
    shrinkageErrorMessage = '';

    try
        [Xos, eta, r_p, k_sav] = optimal_shrinkage_color_fast( ...
            Xhp, opts.shrinkLoss, opts.shrinkKL, opts.shrinkKH); %#ok<ASGLU>

    catch ME
        usedShrinkage = false;
        shrinkageErrorMessage = ME.message;

        warning('optimal_shrinkage_color_fast failed: %s. Falling back to Xhp for distance computation.', ...
            ME.message);

        Xos = Xhp;
        eta = [];
        r_p = [];
        k_sav = [];
    end

    % Pairwise distances on denoised segments
    D = pairwiseEuclidean(Xos);

    % Build artifact template per segment from ORIGINAL segments
    Xhat = zeros(size(Xall));

    for i = 1:nSeg
        d = D(i,:);
        [~, order] = sort(d, 'ascend');
        order(order == i) = []; % exclude self

        Kuse = min(opts.K, numel(order));
        neigh = order(1:Kuse);

        tpl = median(Xall(neigh,:), 1);

        p = segPulseIdx(i);
        win = ones(1, segLen);

        % left overlap
        if p ~= 1
            prevCenter = nearestIndex(timeVals, pulseTimes(p-1));
            thisCenter = nearestIndex(timeVals, pulseTimes(p));
            gap = edPoint - (thisCenter - prevCenter + stPoint);
            if gap > 0 && gap < segLen
                win(1:gap) = sin(pi*(0:gap-1)/(2*gap)).^2;
            end
        end

        % right overlap
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

    % Overlap-add reconstruction
    sa = zeros(nTrials, nSamples);
    z = trialData;

    for i = 1:nSeg
        tr = segTrialIdx(i);
        idx = segSampleIdx(i,:);

        sa(tr, idx) = sa(tr, idx) + Xhat(i,:);
        z(tr, idx) = z(tr, idx) - Xhat(i,:);
    end

    cleanedTrials = z;

    % Standard outputs
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
        assert(~isempty(windowIdx), 'Selected FFT window has no samples.');

        rawWin = trialData(:, windowIdx);
        cleanWin = cleanedTrials(:, windowIdx);

        nWin = numel(windowIdx);
        freqAxis = (0:nWin-1) * (Fs / nWin);

        fftRaw = abs(fft(rawWin, [], 2));
        fftClean = abs(fft(cleanWin, [], 2));

        fftRawMean = mean(fftRaw, 1);
        fftCleanMean = mean(fftClean, 1);
    end

    out = struct();
    out.methodName = sprintf('SMARTAFull_K%d', opts.K);
    out.rawTrials = trialData;
    out.cleanedTrials = cleanedTrials;

    %  Compatibility with demo pipeline
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
    out.diagnostics.usedShrinkage = usedShrinkage;
    out.diagnostics.shrinkageErrorMessage = shrinkageErrorMessage;
    out.diagnostics.hpCutoff = opts.hpCutoff;
    out.diagnostics.shrinkLoss = opts.shrinkLoss;
    out.diagnostics.shrinkKL = opts.shrinkKL;
    out.diagnostics.shrinkKH = opts.shrinkKH;
    out.diagnostics.pulseTimes = pulseTimes;
    out.diagnostics.segmentLength = segLen;
    out.diagnostics.segTrialIdx = segTrialIdx;
    out.diagnostics.segPulseIdx = segPulseIdx;
    out.diagnostics.segSampleIdx = segSampleIdx;
    out.diagnostics.segmentMatrix = Xall;
    out.diagnostics.segmentMatrixHP = Xhp;
    out.diagnostics.segmentMatrixOS = Xos;
    out.diagnostics.artifactSegments = Xhat;
    out.diagnostics.artifactEstimate = sa;
    out.diagnostics.K = opts.K;
    out.diagnostics.eta = eta;
    out.diagnostics.r_p = r_p;
    out.diagnostics.k_sav = k_sav;
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