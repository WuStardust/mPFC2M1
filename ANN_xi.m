close all;clear;clc;
rng('default')
addpath models/
addpath utils/
load("data/data_rat010_0615_spike_train_selected_with_delay.mat")
%% Grid search
ANN_explore_xi = struct( ...
  "H",{}, "xi1",{}, "xi2",{}, "mu",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "s",{}, "Nz",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{} ...
  );
M1Idx = 1; % select M1 neuron
M1spikePart = M1spike(:,M1Idx);
disp('~~~~~~~~~~~~~Start~~~~~~~~~~~~')
tic
parfor i=1:7*7*20
  disp(['===============', num2str(i), ' BEGIN', '===================='])
  s=rng;
  Nz = 5;
  H = 5; % temporal history
  xiList = [10 5 1 0.5 0.1 0.05 0.01];
  xi1 = xiList(ceil(ceil(i/20)/7)); % first stage weight parameters initial range param
  xi2 = xiList(mod(ceil(i/20)-1, 7)+1); % second stage weight parameters initial range param
  mu = 1000; % modified LM algorithm param
  thres = 1e-3; % stop error tolerance
  iterThres = 7; % stop after error over threshold $ times
  maxIter = 500; % max iteration num, over needs re-initial
  alpha = 0; % Regulization parameter
  splitFunc = @(history)splitDataAdvance(1,mPFCspike,M1spikePart,eventTrain,optimalDelay(M1Idx),history);
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
for i=1:7*7*20
  xi1No = ceil(ceil(i/20)/7);
  xi2No = mod(ceil(i/20)-1, 7)+1;
  index = mod(i-1,20)+1;
  ANN_explore_xi(xi1No, xi2No, idx) = results{i};
end
save("results/ANN_explore_xi.mat", "ANN_explore_xi")
