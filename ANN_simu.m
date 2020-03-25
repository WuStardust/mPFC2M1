close all;clear;clc;
rng('default')
addpath models/
addpath utils/
load("data/simu.mat")
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
