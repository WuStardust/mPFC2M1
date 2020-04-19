%% ANN for all neurons
close all;clear;clc;
% rng('default')
addpath models/
addpath utils/
load data/data_rat010_0615_spike_train_top_9.mat
H = 20; % do not explore History
%% ANN
ANN_all_neurons = struct( ...
  "H",{}, "xi1",{}, "xi2",{}, "mu",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "M1No",{}, "s",{}, "Nz",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{}, ...
  "LvalHis",{}, "DBRtrainHis",{}, "DBRvalHis",{}, "Whis",{},"muHis",{},"HessianDetHis",{} ...
  );
M1count = 9;
repeatNum = 20;
disp('~~~~~~~~~~~~~Start~~~~~~~~~~~~')
tic
for M1Idx=1:M1count
  disp(['++++++M1 Idx: ', num2str(M1Idx), '++++++'])
  % extract params
  M1No = M1index(M1Idx); % M1 neuron No.
  mPFCspikePart = mPFCspike(:,mPFCm1Map(M1Idx,:)); % corresponding mPFC neurons
  M1spikePart = M1spike(:,M1Idx); % current M1 neuron
  delay = M1delay(M1Idx);
  for Nz=5:15
  parfor i=1:repeatNum
    disp(['===============', datestr(datetime), '-', num2str(Nz), '-', num2str(i), '==============='])
    s=rng;
    xi1 = 0.1; % first stage weight parameters initial range param
    xi2 = 0.5; % second stage weight parameters initial range param
    mu = 1000; % modified LM algorithm param
    thres = 1e-3; % stop error tolerance
    iterThres = 7; % stop after error over threshold $ times
    maxIter = 1000; % max iteration num, over needs re-initial
    alpha = 0; % Regulization parameter
    splitFunc = @(history)splitDataAdvance(1,mPFCspikePart,M1spikePart,eventTrain,delay,segTrain,history);
    verbose = 2;
    [W,L,DBR,Lval,LtrainHis,LvalHis,DBRtrainHis,DBRvalHis,Whis,muHis,HessianDetHis] = runANN(H, Nz, xi1, xi2, mu, thres, iterThres, maxIter, alpha, splitFunc, verbose);
    ANN_all_neurons(M1Idx, Nz, i) = struct( ...
      "H",H, "xi1",xi1, "xi2",xi2, "mu",mu, "thres",thres, "iterThres",iterThres, ...
      "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, "M1No",M1No, "s",s, "Nz",Nz, ...
      "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LtrainHis, ...
      "LvalHis",LvalHis, "DBRtrainHis",DBRtrainHis, "DBRvalHis",DBRvalHis, "Whis",Whis, "muHis",muHis,"HessianDetHis",HessianDetHis ...
      );
  end
  end
  ANN_neuron = squeeze(ANN_all_neurons(M1Idx,:,:));
  save(['results/ANN_top_9_neurons/', num2str(M1Idx), '.mat'], "ANN_neuron")
end
toc
disp('~~~~~~~~~~~~~End~~~~~~~~~~~~')
save("results/ANN_top_9_neurons.mat", "ANN_all_neurons")
