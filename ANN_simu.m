close all;clear;clc;
rng('default')
addpath models/
addpath utils/
load("data/simu2nd.mat")
input=input';
output=output';
%% Grid search
disp(['===============', datestr(datetime), '==============='])
Nz = 15;
H = 50; % temporal history, todo: grid search
xi1 = 0.1; % first stage weight parameters initial range param
xi2 = 0.5; % second stage weight parameters initial range param
mu = 1; % modified LM algorithm param
thres = 1e-3; % stop error tolerance
iterThres = 7; % stop after error over threshold $ times
maxIter = 1000; % max iteration num, over needs re-initial
alpha = 0; % Regulization parameter
splitFunc = @(history)splitData(input,output,history);
verbose = 0;
[W, L, DBR, Lval, LHistory] = runANN(H, Nz, xi1, xi2, mu, thres, iterThres, maxIter, alpha, splitFunc, verbose);
% for i=1:800
%   [Nz,idx] = getParamIndex(i);
%   ANN_explore_Nz(Nz, idx) = results{i};
% end
% save("results/ANN_explore_Nz.mat", "ANN_explore_Nz")
%%
H = 50; % temporal history
xi=0.1; % Weights random initial range parameter
alpha=0; % Regulization parameter
thres = 1e-3; % stop error tolerance
iterThres = 0; % stop after error over threshold $ times
maxIter = 40; % max iteration num, over needs re-initial
M1Idx = 2; % select M1 neuron
splitFun = @(history)splitDataLocal(input,output,history, 5, 2, 0.5);
verbose = 0;
[W, L, GLM_sec_DBR] = runGLM(H, xi, thres, iterThres, maxIter, alpha, splitFun, verbose);

%%
function [trainLen, valLen, testLen, trainX, valX, testX, trainY, valY, testY] = splitDataLocal(inputSignal,outputSignal,H, lagnum, lagOrder, lagalpha)
%Split data into train, val, test
Xhat = ensembleLaguerre(inputSignal, H, lagnum, lagOrder, lagalpha);
trainLen = ceil(0.6*length(inputSignal));
valLen   = ceil(0.2*length(inputSignal));
testLen  = length(Xhat) - trainLen - valLen;

trainX = Xhat(1:trainLen,:);
valX   = Xhat(trainLen+1:trainLen+valLen,:);
testX  = Xhat(trainLen+valLen+1:end,:);
trainY = outputSignal(H:H+trainLen-1);
valY   = outputSignal(H+trainLen:H+trainLen+valLen-1);
testY  = outputSignal(H+trainLen+valLen:end);
end
