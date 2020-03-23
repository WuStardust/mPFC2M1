close all;clear;clc;
rng('default')
addpath models/
addpath utils/
load("data/data_rat010_0615_spike_train_selected_with_delay.mat")
%% Grid search
ANN_explore_Nz = struct( ...
  "H",{}, "xi1",{}, "xi2",{}, "mu",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "s",{}, "Nz",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{} ...
  );
M1Idx = 1; % select M1 neuron
M1spikePart = M1spike(:,M1Idx);
disp('~~~~~~~~~~~~~Start~~~~~~~~~~~~')
tic
for Nz=1:25
parfor i=1:40
  disp(['===============', num2str(i), ' BEGIN', '===================='])
  s=rng;
%   Nz = getParamIndex(i);
  H = 15; % temporal history, todo: grid search
  xi1 = 0.1; % first stage weight parameters initial range param
  xi2 = 0.5; % second stage weight parameters initial range param
  mu = 1000; % modified LM algorithm param
  thres = 1e-3; % stop error tolerance
  iterThres = 7; % stop after error over threshold $ times
  maxIter = 1000; % max iteration num, over needs re-initial
  alpha = 0; % Regulization parameter
  splitFunc = @(history)splitDataAdvance(1,mPFCspike,M1spikePart,eventTrain,optimalDelay(M1Idx),segTrain,history);
  verbose = 2;
  [W, L, DBR, Lval, LHistory] = runANN(H, Nz, xi1, xi2, mu, thres, iterThres, maxIter, alpha, splitFunc, verbose);
  ANN_explore_Nz(Nz, i) = struct( ...
    "H",H, "xi1",xi1, "xi2",xi2, "mu",mu, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, "s",s, "Nz",Nz, ...
    "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LHistory ...
    );
  if (i/10 == 0)
    disp(['      Nz ', num2str(Nz, '%02d'), ' idx ', num2str(mod(i-1, repeatNum)+1, '%02d'), ' Lval ', num2str(Lval)])
  end
end
end
toc
% for i=1:800
%   [Nz,idx] = getParamIndex(i);
%   ANN_explore_Nz(Nz, idx) = results{i};
% end
save("results/ANN_explore_Nz.mat", "ANN_explore_Nz")
