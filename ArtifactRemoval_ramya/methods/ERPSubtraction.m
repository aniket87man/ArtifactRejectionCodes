function out = ERPSubtraction(trialData,timeVals,params)
% ERPSubtraction
%
% Baseline artifact-removal method
% compute ERP/template across selected trials and subtract it from each trial.
%
% Inputs
%   trialData : trials x time
%   timeVals  : 1 x time, in seconds
%   params    : optional struct
%
% Optional params
%   params.subtractWindow        : [t1 t2], default = full trial
%   params.doBaselineCorrection  : true/false, default = false
%   params.baselineWindow        : [t1 t2], default = [-0.7 -0.2]
%
% Output
%   out.cleanedData
%   out.template
%   out.fullERP
%   out.rawMean
%   out.cleanedMean
%   out.methodName
%   out.params

if nargin < 3
    params = struct();
end

if ~isfield(params,'subtractWindow') || isempty(params.subtractWindow)
    params.subtractWindow = [timeVals(1) timeVals(end)];
end

if ~isfield(params,'doBaselineCorrection') || isempty(params.doBaselineCorrection)
    params.doBaselineCorrection = false;
end

if ~isfield(params,'baselineWindow') || isempty(params.baselineWindow)
    params.baselineWindow = [-0.7 -0.2];
end

% Safety checks
if size(trialData,2) ~= length(timeVals)
    error('trialData must be trials x time, and length(timeVals) must match size(trialData,2).');
end

dataIn = trialData;

% Optional baseline correction
if params.doBaselineCorrection
    baselineIdx = timeVals > params.baselineWindow(1) & timeVals < params.baselineWindow(2);
    if ~any(baselineIdx)
        error('No samples found in baselineWindow.');
    end
    baselineVals = mean(dataIn(:,baselineIdx),2);
    dataIn = dataIn - baselineVals;
end

% ERP/template
fullERP = mean(dataIn,1);

% Subtract either full ERP or windowed ERP
subtractIdx = timeVals > params.subtractWindow(1) & timeVals <= params.subtractWindow(2);
if ~any(subtractIdx)
    error('No samples found in subtractWindow.');
end

template = zeros(size(fullERP));
template(subtractIdx) = fullERP(subtractIdx);

cleanedData = dataIn - template;

% Output
out = struct();
out.methodName = 'ERP subtraction baseline';
out.cleanedData = cleanedData;
out.template = template;
out.fullERP = fullERP;
out.rawMean = mean(dataIn,1);
out.cleanedMean = mean(cleanedData,1);
out.timeVals = timeVals;
out.params = params;
end