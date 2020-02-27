%% ANN

close all;clear;clc;
rng('default')
addpath models/
addpath utils/
load("data/data_rat010_0615_spike_train_selected.mat")
%% Grid search
ANN_explore_H = struct( ...
  "H",{}, "xi1",{}, "xi2",{}, "mu",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "s",{}, "Nz",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{} ...
  );
HSearchNum = 10;
repeatNum = 100;
M1Idx = 1; % select M1 neuron
M1spikePart = M1spike(:,M1Idx);
tic
parfor i=1:HSearchNum*repeatNum
  Hlist = 1:HSearchNum;
  s=rng;
  H = Hlist(ceil(i/repeatNum));
  Nz = 10; % temporal history
  xi1 = 0.05; % first stage weight parameters initial range param
  xi2 = 0.1; % second stage weight parameters initial range param
  mu = 1000; % modified LM algorithm param
  thres = 1e-3; % stop error tolerance
  iterThres = 7; % stop after error over threshold $ times
  maxIter = 1000; % max iteration num, over needs re-initial
  alpha = 0; % Regulization parameter
  splitFunc = @(history)splitData(mPFCspike,M1spikePart,history); % choose splitData function
  verbose = 2;
  [W, L, DBR, Lval, LHistory] = runANN(H, Nz, xi1, xi2, mu, thres, iterThres, maxIter, alpha, splitFunc, verbose);
  results{i} = struct( ...
    "H",H, "xi1",xi1, "xi2",xi2, "mu",mu, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, "s",s, "Nz",Nz, ...
    "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LHistory ...
    );
  if (i/10 == 0)
    disp(['      Nz ', num2str(Nz, '%02d'), ' idx ', num2str(mod(i-1, repeatNum)+1, '%02d'), ' Lval ', num2str(Lval)])
  end
end
toc
Hlist = 1:HSearchNum;
for i=1:HSearchNum*repeatNum
  H = Hlist(ceil(i/repeatNum));
  idx = mod(i-1, repeatNum)+1;
  ANN_explore_H(H, idx) = results{i};
end
save("results/ANN_explore_H.mat", "ANN_explore_H")
