%% ANN

close all;clear;clc;
rng('default')
addpath models/
addpath utils/
load("data/data_rat010_0615_spike_train_selected.mat")
%% Grid search
% explore Nz as H=5, rng use default
% 
% for Nz=8, 10, He is singular, need re-initial. 
% 
% result seems not good, pretrain needed?

ANN_explore_Nz = struct( ...
  "H",{}, "xi1",{}, "xi2",{}, "mu",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "s",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{} ...
  );
NzSearchNum = 20;
repeatNum = 20;
M1Idx = 1; % select M1 neuron
M1spikePart = M1spike(:,M1Idx);
tic
parfor i=1:NzSearchNum*repeatNum
  Nzlist = 1:20;
  s=rng;
  Nz = Nzlist(ceil(i/20));
  H = 5; % temporal history, todo: grid search
  xi1 = 0.05; % first stage weight parameters initial range param
  xi2 = 0.1; % second stage weight parameters initial range param
  mu = 1000; % modified LM algorithm param
  thres = 1e-3; % stop error tolerance
  iterThres = 7; % stop after error over threshold $ times
  maxIter = 5000; % max iteration num, over needs re-initial
  alpha = 0; % Regulization parameter
  splitFunc = @(history)splitData(mPFCspike,M1spikePart,history); % choose splitData function
  verbose = 3;
  [W, L, DBR, Lval, LHistory] = runANN(H, Nz, xi1, xi2, mu, thres, iterThres, maxIter, alpha, splitFunc, verbose);
  results{i} = struct( ...
    "H",H, "xi1",xi1, "xi2",xi2, "mu",mu, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, "s",s, ...
    "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LHistory ...
    );
  if (i/10 == 0)
    disp(['      Nz ', num2str(Nz, '%02d'), ' idx ', num2str(mod(i-1, 20)+1, '%02d'), ' Lval ', num2str(Lval)])
  end
end
toc
Nzlist = 1:20;
for i=1:NzSearchNum*repeatNum
  Nz = Nzlist(ceil(i/20));
  idx = mod(i-1, 20)+1;
  ANN_explore_Nz(Nz, idx) = results{i};
end
save("results/ANN_explore_Nz.mat", "ANN_explore_Nz")