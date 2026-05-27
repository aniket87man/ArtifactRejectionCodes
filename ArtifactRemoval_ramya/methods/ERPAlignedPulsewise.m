function out = ERPAlignedPulsewise(trialData,timeVals,params)
% ERPAlignedPulsewise
%
% Hybrid artifact-removal method.
%
% Step 1:
%   Apply aligned ERP template subtraction to remove the broad
%   stimulation-locked artifact component.
%
% Step 2:
%   Apply pulse-wise template subtraction to the ERPAligned residual to
%   remove remaining pulse-locked residual structure.
%
% Inputs
%   trialData : trials x time
%   timeVals  : 1 x data points,time, seconds
%   params    : struct
%
% Optional params
%   params.subtractWindow : [t1 t2], default = [0 0.4]
%   params.alignedParams  : params passed to ERPAligned
%   params.pulseParams    : params passed to PulsewiseTemplate
%
% Output
%   out.cleanedData
%   out.erpAlignedOut
%   out.pulseOut
%   out.totalTemplatePerTrial
%   out.rawMean
%   out.afterERPAlignedMean
%   out.cleanedMean
%   out.params

if nargin < 3
    params = struct();
end

if size(trialData,2) ~= length(timeVals)
    error('trialData must be trials x time, and length(timeVals) must match size(trialData,2).');
end

%% Common subtraction window
if ~isfield(params,'subtractWindow') || isempty(params.subtractWindow)
    params.subtractWindow = [0 0.4];
end

%% ERPAligned parameters
if ~isfield(params,'alignedParams') || isempty(params.alignedParams)

    alignedParams = struct();
    alignedParams.subtractWindow = params.subtractWindow;
    alignedParams.alignWindow = [-0.01 0.03];
    alignedParams.maxShiftMS = 10;
    alignedParams.doBaselineCorrection = false;

else
    alignedParams = params.alignedParams;

    if ~isfield(alignedParams,'subtractWindow') || isempty(alignedParams.subtractWindow)
        alignedParams.subtractWindow = params.subtractWindow;
    end
end

%% PulsewiseTemplate parameters
if ~isfield(params,'pulseParams') || isempty(params.pulseParams)

    pulseParams = struct();

    pulseParams.subtractWindow = params.subtractWindow;

    pulseParams.pulseSearchWindow = [-0.02 0.32];
    pulseParams.pulseWindowMS = [-5 35];
    pulseParams.minPulseDistanceMS = 25;
    pulseParams.expectedNumPulses = 7;
    pulseParams.maxNumPulses = 20;

    % Lower threshold is used because this method detects pulse residuals
    % after ERPAligned subtraction.
    pulseParams.thresholdMAD = 3;

    pulseParams.localBaselineMS = 20;
    pulseParams.templateStatistic = 'median';
    pulseParams.taperEdgeMS = 2;
    pulseParams.doBaselineCorrection = false;

else
    pulseParams = params.pulseParams;

    if ~isfield(pulseParams,'subtractWindow') || isempty(pulseParams.subtractWindow)
        pulseParams.subtractWindow = params.subtractWindow;
    end
end

%% Step 1: aligned ERP template subtraction
erpAlignedOut = ERPAligned(trialData,timeVals,alignedParams);

residualAfterERPAligned = erpAlignedOut.cleanedData;

%% Step 2: pulse-wise template subtraction on ERPAligned residual
pulseOut = PulsewiseTemplate(residualAfterERPAligned,timeVals,pulseParams);

cleanedData = pulseOut.cleanedData;

%% Total template
% ERPAligned uses a trial-specific template.
% PulsewiseTemplate uses one residual pulse-wise template for all trials.
numTrials = size(trialData,1);

pulseTemplatePerTrial = repmat(pulseOut.fullTemplate,numTrials,1);
totalTemplatePerTrial = erpAlignedOut.trialTemplates + pulseTemplatePerTrial;

%% Output
out = struct();

out.methodName = 'ERPAligned plus PulsewiseTemplate';

out.cleanedData = cleanedData;

out.erpAlignedOut = erpAlignedOut;
out.pulseOut = pulseOut;

out.totalTemplatePerTrial = totalTemplatePerTrial;
out.pulseTemplatePerTrial = pulseTemplatePerTrial;

out.rawMean = mean(trialData,1);
out.afterERPAlignedMean = mean(residualAfterERPAligned,1);
out.cleanedMean = mean(cleanedData,1);

out.timeVals = timeVals;
out.params = params;
out.alignedParams = alignedParams;
out.pulseParams = pulseParams;

end