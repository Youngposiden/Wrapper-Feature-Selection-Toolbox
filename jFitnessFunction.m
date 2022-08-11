% Fitness Function KNN (9/12/2020)

function cost = jFitnessFunction(feat,label,X,opts)
% Default of [alpha; beta]
ws = [0.99; 0.01];

if isfield(opts,'ws'), ws = opts.ws; end

% Check if any feature exist
if sum(X == 1) == 0
  cost = 1;
else
  % Error rate
    switch opts.method % will be slow, should choose on initalization and eliminate
        case 'knn'
            error    = jwrapper_KNN(feat(:,X == 1),label,opts);
        case 'dt'
            error = jwrapper_DT(feat(:,X == 1),label,opts);
        case 'rf'
            error = jwrapper_RF(feat(:,X == 1),label,opts);
    end
  % Number of selected features
  num_feat = sum(X == 1);
  % Total number of features
  max_feat = length(X); 
  % Set alpha & beta
  alpha    = ws(1); 
  beta     = ws(2);
  % Cost function 
  cost     = alpha * error + beta * (num_feat / max_feat); 
end
end


%% ---Call Functions-----------------------------------------------------
% may want to move into seperate files for ease of result identification,
% unless the local functions have a significant impact on speed 
function error = jwrapper_KNN(sFeat,label,opts)
if isfield(opts,'k'), k = opts.k; end
if isfield(opts,'Model'), Model = opts.Model; end

% Define training & validation sets
trainIdx = Model.training;    testIdx = Model.test;
xtrain   = sFeat(trainIdx,:); ytrain  = label(trainIdx);
xvalid   = sFeat(testIdx,:);  yvalid  = label(testIdx);
% Training model
My_Model = fitcknn(xtrain,ytrain,'NumNeighbors',k); 
% Prediction
pred     = predict(My_Model,xvalid);
% Accuracy
Acc      = sum(pred == yvalid) / length(yvalid);
% Error rate
error    = 1 - Acc; 
end

function error = jwrapper_DT(sFeat,label,opts)
if isfield(opts,'Model'), Model = opts.Model; end

% Define training & validation sets
trainIdx = Model.training;    testIdx = Model.test;
xtrain   = sFeat(trainIdx,:); ytrain  = label(trainIdx);
xvalid   = sFeat(testIdx,:);  yvalid  = label(testIdx);
% Training model
My_Model = fitctree(xtrain,ytrain); 
% Prediction
pred     = predict(My_Model,xvalid);
% Accuracy
Acc      = sum(pred == yvalid) / length(yvalid);
% Error rate
error    = 1 - Acc; 
end

%% slow and painful (training a lot of trees)
function error = jwrapper_RF(sFeat,label,opts)
if isfield(opts,'Model'), Model = opts.Model; end
if isfield(opts,'numtrees'), numtrees = opts.numtrees; end

% Define training & validation sets
trainIdx = Model.training;    testIdx = Model.test;
xtrain   = sFeat(trainIdx,:); ytrain  = label(trainIdx);
xvalid   = sFeat(testIdx,:);  yvalid  = label(testIdx);

% Training model
t = templateTree('NumVariablesToSample','all','PredictorSelection','interaction-curvature','Surrogate','on');
My_Model = fitcensemble(xtrain,ytrain,'Method','Bag','NumLearningCycles',numtrees, 'Learners',t);
% Prediction
pred     = predict(My_Model,xvalid);
% Accuracy
Acc      = sum(pred == yvalid) / length(yvalid);
% Error rate
error    = 1 - Acc; 
end