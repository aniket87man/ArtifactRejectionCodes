function out = PCATemplate(trialData,timeVals,params)
% PCATemplate
%
% PCA-based stimulation artifact template subtraction.
%
% Main idea:
% 1. Select an artifact window.
% 2. Compute the mean artifact template in that window.
% 3. Compute PCA/SVD on trial-to-trial residuals in that window.
% 4. Reconstruct the top K artifact components.
% 5. Subtract mean template + top K PCA components from each trial.
%
% Inputs
% trialData : trials x time
% timeVals  : 1 x time, seconds
% params    : struct
%
% Optional params
% params.artifactWindow        : [t1 t2], default = [0 0.4]
% params.numComponents         : number of PCA components, default = 3
% params.removeMeanTemplate    : true/false, default = true
% params.taperEdgeMS           : taper edge of artifact window, default = 2 ms
% params.doBaselineCorrection  : true/false, default = false
% params.baselineWindow        : [t1 t2], default = [-0.7 -0.2]
%
% Output
% out.cleanedData
% out.artifactModel
% out.meanTemplate
% out.pcTemplates
% out.scores
% out.explainedVariance
% out.artifactIdx
% out.rawMean
% out.cleanedMean
% out.params

if nargin < 3
    params = struct();
end

if ~isfield(params,'artifactWindow') || isempty(params.artifactWindow)
    params.artifactWindow = [0 0.4];
end

if ~isfield(params,'numComponents') || isempty(params.numComponents)
    params.numComponents = 3;
end

if ~isfield(params,'removeMeanTemplate') || isempty(params.removeMeanTemplate)
    params.removeMeanTemplate = true;
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
    baselineIdx = timeVals >= params.baselineWindow(1) & timeVals < params.baselineWindow(2);

    if ~any(baselineIdx)
        error('No samples found in baselineWindow.');
    end

    baselineVals = mean(dataIn(:,baselineIdx),2);
    dataIn = dataIn - baselineVals;
end

Fs = 1/median(diff(timeVals));

artifactIdx = find(timeVals > params.artifactWindow(1) & ...
                   timeVals <=  params.artifactWindow(2));

if isempty(artifactIdx)
    error('No samples found in artifactWindow.');
end

numTrials = size(dataIn,1);
numTime = size(dataIn,2);

%% Extract artifact-window data
X = dataIn(:,artifactIdx);

%% Mean artifact template
meanTemplateWindow = mean(X,1);

if params.removeMeanTemplate
    Xcentered = X - meanTemplateWindow;
else
    Xcentered = X;
    meanTemplateWindow = zeros(size(meanTemplateWindow));
end

%% PCA using SVD
% Xcentered = U*S*V'
[U,S,V] = svd(Xcentered,'econ');

singularValues = diag(S);
componentPower = singularValues.^2;
explainedVariance = componentPower ./ (sum(componentPower) + eps);

maxComponents = min([params.numComponents,size(V,2),size(U,2)]);

if maxComponents < 1
    pcaModelCentered = zeros(size(X));
    pcTemplates = [];
    scores = [];
else
    scores = U(:,1:maxComponents) * S(1:maxComponents,1:maxComponents);
    pcTemplates = V(:,1:maxComponents)';

    pcaModelCentered = scores * V(:,1:maxComponents)';
end

%% Artifact model in artifact window
artifactModelWindow = meanTemplateWindow + pcaModelCentered;

%% Taper artifact model edges to reduce discontinuities
edgeSamples = round((params.taperEdgeMS/1000) * Fs);
artifactModelWindow = applyWindowTaper(artifactModelWindow,edgeSamples);

%% Full artifact model
artifactModel = zeros(numTrials,numTime);
artifactModel(:,artifactIdx) = artifactModelWindow;

%% Subtract artifact model
cleanedData = dataIn - artifactModel;

%% Output
out = struct();

out.methodName = 'PCA template subtraction';

out.cleanedData = cleanedData;
out.artifactModel = artifactModel;

out.meanTemplate = zeros(1,numTime);
out.meanTemplate(artifactIdx) = meanTemplateWindow;

out.pcTemplates = pcTemplates;
out.scores = scores;
out.explainedVariance = explainedVariance;

out.artifactIdx = artifactIdx;
out.artifactWindow = params.artifactWindow;

out.rawMean = mean(dataIn,1);
out.cleanedMean = mean(cleanedData,1);

out.timeVals = timeVals;
out.params = params;

end

function Y = applyWindowTaper(X,edgeSamples)
% Apply raised-cosine taper to the edges of each row of X.

Y = X;

if edgeSamples <= 0
    return;
end

n = size(X,2);
edgeSamples = min(edgeSamples,floor(n/2));

if edgeSamples <= 0
    return;
end

t = linspace(0,pi/2,edgeSamples);
leftTaper = sin(t).^2;
rightTaper = fliplr(leftTaper);

Y(:,1:edgeSamples) = Y(:,1:edgeSamples) .* leftTaper;
Y(:,end-edgeSamples+1:end) = Y(:,end-edgeSamples+1:end) .* rightTaper;

end