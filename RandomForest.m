%% YQZ, 2020.07.27
% Generating a dummy random forest function for MPI test

function [ ResponseV,Mdl] = RandomForest(INPUTS)

%% Number of predictor variables
% PredictorNumber = 10;

%% generate dummy data


X         = INPUTS(:,2:size(INPUTS,2));
Y         = INPUTS(:,1);

%% Establish regression tree ensemble (bagging, a generalised random forest method)
t         = templateTree('NumPredictorsToSample','all',...
          'PredictorSelection','interaction-curvature','Surrogate','on');
Mdl       = fitrensemble(X,Y,'Method','Bag','NumLearningCycles',100,...
           'Learners',t);
       
% paroptions = statset('UseParallel',true);

% Mdl       = fitrensemble(X,Y,'Method','Bag',...
%            'Learners',t,'OptimizeHyperparameters','auto');
       
%% optimise the OptimizeHyperparameters within Random forest
% t   = templateTree('Surrogate','on');
% Mdl = fitrensemble(X,Y,'Learners',t,...
%     'OptimizeHyperparameters',{'NumLearningCycles','LearnRate','MaxNumSplits'},'HyperparameterOptimizationOptions',struct('UseParallel',true));%,'UseParallel','true'
%% Prediction and cross validations
yPre      = predict(Mdl,Mdl.X);
CVMdl     = crossval(Mdl,'kfold',10);
yCrossVal = kfoldPredict(CVMdl);
%% Outputing observed and predicted response variable

ResponseV = [Mdl.Y yPre yCrossVal]; 

end

