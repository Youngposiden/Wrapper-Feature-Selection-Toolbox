function Acc = jrf(sFeat,label,opts)
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

fprintf('\n Accuracy: %g %%',100 * Acc);
end