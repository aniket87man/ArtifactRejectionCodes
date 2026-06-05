function fftOut = compute_fft_summary(trialData,timeVals,windowIdx)
% compute_fft_summary
%
% Computes FFT magnitude for trials x time data.
%
% Inputs
%   trialData 
%   timeVals 
%   windowIdx 
%
% Output
%   fftOut.freqAxis
%   fftOut.meanMagnitude
%   fftOut.logMeanMagnitude
%   fftOut.allMagnitude
%   fftOut.Fs

if size(trialData,2) ~= length(timeVals)
    error('trialData must be trials x time, and length(timeVals) must match size(trialData,2).');
end

Fs = 1/median(diff(timeVals));

dataWindow = trialData(:,windowIdx);
nFFT = size(dataWindow,2);

X = fft(dataWindow,[],2);
mag = abs(X);

freqAxis = (0:nFFT-1) * Fs/nFFT;

fftOut = struct();
fftOut.freqAxis = freqAxis;
fftOut.allMagnitude = mag;
fftOut.meanMagnitude = mean(mag,1);
fftOut.logMeanMagnitude = log10(fftOut.meanMagnitude + eps);
fftOut.Fs = Fs;
fftOut.nFFT = nFFT;
end