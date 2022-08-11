function Acc = jdt(sFeat,label,opts)
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

fprintf('\n Accuracy: %g %%',100 * Acc);
end