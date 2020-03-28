%% ANN

close all;clear;clc;
rng('default')
addpath models/
addpath utils/
load("data/data_rat010_0615_spike_train_selected_with_delay.mat")
%% Grid search
ANN_explore_H = struct( ...
  "H",{}, "xi1",{}, "xi2",{}, "mu",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "s",{}, "Nz",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{}, ...
  "LvalHis",{}, "DBRtrainHis",{}, "DBRvalHis",{}, "Whis",{} ...
  );
HSearchNum = 20;
repeatNum = 30;
M1Idx = 3; % select M1 neuron
M1spikePart = M1spike(:,M1Idx);
disp('~~~~~~~~~~~~~Start~~~~~~~~~~~~')
tic
for H=1:1:HSearchNum
parfor i=1:repeatNum
  disp(['===============', datestr(datetime), '-', num2str(H), '-', num2str(i), '==============='])
  s=rng;
  Nz = 10; % hidden neuron number
  xi1 = 0.1; % first stage weight parameters initial range param
  xi2 = 0.5; % second stage weight parameters initial range param
  mu = 1000; % modified LM algorithm param
  thres = 1e-3; % stop error tolerance
  iterThres = 7; % stop after error over threshold $ times
  maxIter = 1000; % max iteration num, over needs re-initial
  alpha = 0; % Regulization parameter
  splitFunc = @(history)splitDataAdvance(1,mPFCspike,M1spikePart,eventTrain,optimalDelay(M1Idx),segTrain,history);
  verbose = 2;
  [W,L,DBR,Lval,LtrainHis,LvalHis,DBRtrainHis,DBRvalHis,Whis] = runANN(H, Nz, xi1, xi2, mu, thres, iterThres, maxIter, alpha, splitFunc, verbose);
  ANN_explore_H(H, i) = struct( ...
    "H",H, "xi1",xi1, "xi2",xi2, "mu",mu, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, "s",s, "Nz",Nz, ...
    "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LtrainHis, ...
    "LvalHis",LvalHis, "DBRtrainHis",DBRtrainHis, "DBRvalHis",DBRvalHis, "Whis",Whis ...
    );
end
end
toc
disp('~~~~~~~~~~~~~End~~~~~~~~~~~~')
% Hlist = 1:HSearchNum;
% for i=1:HSearchNum*repeatNum
%   H = Hlist(ceil(i/repeatNum));
%   idx = mod(i-1, repeatNum)+1;
%   ANN_explore_H(H, idx) = results{i};
% end
save("results/ANN_explore_H_3.mat", "ANN_explore_H")
