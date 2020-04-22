close all;clear;clc;
rng('default')
addpath models\
addpath utils\
load data\data_rat010_0615_spike_train_top_9_mPFC_8.mat
H = 20; % do not explore History
%% GLM
GLM_all_neurons = struct( ...
  "H",{}, "xi",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "M1No",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{} ...
  );

for M1Idx=1:9 % For each M1 neuron
  disp(['=====M1 Idx: ', num2str(M1Idx), '====='])
  % extract params
  M1No = M1index(M1Idx); % M1 neuron No.
  mPFCspikePart = mPFCspike(:,mPFCm1Map(M1Idx,:)); % corresponding mPFC neurons
  M1spikePart = M1spike(:,M1Idx); % current M1 neuron
  delay = M1delay(M1Idx);
  % hyperparams
  xi=0.1; % Weights random initial range parameter
  alpha=0; % Regulization parameter
  thres = 1e-3; % stop error tolerance
  iterThres = 0; % stop after error over threshold $ times
  maxIter = 40; % max iteration num, over needs re-initial
  % run the model
  splitFun = @(history)splitDataAdvance(1,mPFCspikePart,M1spikePart,eventTrain,delay,segTrain,history); % choose splitData function as first Order
  verbose = 2;
  [W, L, DBR, Lval, LHistory] = runGLM(H, xi, thres, iterThres, maxIter, alpha, splitFun, verbose);
  GLM_all_neurons(M1Idx) = struct( ...
    "H",H, "xi",xi, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, "M1No",M1No, ...
    "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LHistory ...
    );
end
%% GLM 2nd order
GLM_sec_all_neurons = struct( ...
  "H",{}, "xi",{}, "thres",{}, "iterThres",{}, "lagnum", {}, "lagalpha", {}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, "M1No",{}, ...
  "W",{}, "L",{}, "DBR",{} ...
  );
for M1Idx=1:9 % For each M1 neuron
  disp(['=====M1 Idx: ', num2str(M1Idx), '====='])
  % extract params
  M1No = M1index(M1Idx); % M1 neuron No.
  mPFCspikePart = mPFCspike(:,mPFCm1Map(M1Idx,:)); % corresponding mPFC neurons
  M1spikePart = M1spike(:,M1Idx); % current M1 neuron
  delay = M1delay(M1Idx);
  bestDBR=10;
  for lagnum=1:min(H-1,5)
    for lagalpha=0.05:0.05:0.7
      xi=0.1; % Weights random initial range parameter
      alpha=0; % Regulization parameter
      thres = 1e-3; % stop error tolerance
      iterThres = 0; % stop after error over threshold $ times
      maxIter = 40; % max iteration num, over needs re-initial
      order = ['vm2-',num2str(lagnum),'-',num2str(lagalpha)];
      splitFun = @(history)splitDataAdvance(order,mPFCspikePart,M1spikePart,eventTrain,delay,segTrain,history);
      verbose = 2;
      disp(['order: ', order])
      [W, L, GLM_sec_DBR] = runGLM(H, xi, thres, iterThres, maxIter, alpha, splitFun, verbose);
      if (GLM_sec_DBR<bestDBR)
        bestW = W;
        bestL = L;
        bestDBR = GLM_sec_DBR;
        bestlagnum = lagnum;
        bestlagalpha = lagalpha;
      end
    end
  end
  GLM_sec_all_neurons(M1Idx) = struct( ...
    "H",H, "xi",xi, "thres",thres, "iterThres",iterThres, ...
    "lagnum", bestlagnum, "lagalpha", bestlagalpha, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, "M1No",M1No, ...
    "W",bestW, "L",bestL, "DBR",bestDBR ...
    );
end
%% save results
save("results/baseline-top9-neurons-mPFC-8.mat", "GLM_all_neurons", "GLM_sec_all_neurons");
%% get statistical
% load data; uncoment following lines when needed
close all;clear;clc;
addpath models/
load data/data_rat010_0615_spike_train_top_9_mPFC_8.mat
load results/baseline-top9-neurons-mPFC-8.mat
neurons_DBR = zeros(3, 9);
for i=1:9
  neurons_DBR(1, i) = GLM_all_neurons(i).DBR;
  neurons_DBR(2, i) = GLM_sec_all_neurons(i).DBR;
end

neurons_cc = zeros(3, 9);
for i=1:9
  mPFCspikePart = mPFCspike(:,mPFCm1Map(i,:)); % corresponding mPFC neurons
  M1spikePart = M1spike(:,i); % current M1 neuron
  delay = M1delay(i);
  % GLM
  W     = GLM_all_neurons(i).W;
  H     = GLM_all_neurons(i).H;
  [~,~,~,~,~,testX,~,~,~] = splitDataAdvance(1,mPFCspikePart,M1spikePart,eventTrain,delay,segTrain,H);
  GLMtestLambdaYpre = GLMmodel(testX, W);
  % GLM sec
  W        = GLM_sec_all_neurons(i).W;
  H        = GLM_sec_all_neurons(i).H;
  lagnum   = GLM_sec_all_neurons(i).lagnum;
  lagalpha = GLM_sec_all_neurons(i).lagalpha;
  order = ['vm2-',num2str(lagnum),'-',num2str(lagalpha)];
  [~,~,~,~,~,testX,~,~,testY] = splitDataAdvance(order,mPFCspikePart,M1spikePart,eventTrain,delay,segTrain,H);
  GLMsectestLambdaYpre = GLMmodel(testX, W);
  
  kernelSize = 100;
  smoothedLambda = gaussianSmooth(testY, kernelSize);
  smoothedGLMLambdaPre = gaussianSmooth(GLMtestLambdaYpre, kernelSize);
  smoothedGLMsecLambdaPre = gaussianSmooth(GLMsectestLambdaYpre, kernelSize);
  
  cc = corrcoef(smoothedLambda, smoothedGLMLambdaPre);
  neurons_cc(1,i) = cc(2);
  cc = corrcoef(smoothedLambda, smoothedGLMsecLambdaPre);
  neurons_cc(2,i) = cc(2);

%   % plot
%   t = 200:0.01:300;
%   index = 20000:30000;
%   h = figure("Name", ['top', num2str(i), 'neuron']);
%   hold on
%   plot(t, smoothedLambda(index), 'r');
%   plot(t, smoothedGLMLambdaPre(index), 'b');
%   plot(t, smoothedGLMsecLambdaPre(index), 'g');
%   hold off
%   legend("Actual M1 spike", "GLM", "2nd-Order GLM", "Box","off", "Orientation","horizontal")
%   savefig(h, ['results/top9/top', num2str(i), 'neuron.fig'])
end
%% save statistical results
h = figure("Name", "Neurons-DBR-cc");
subplot(1,2,1)
hold on
bar(neurons_DBR');
plot(0.5:0.01:9.3, ones(1, 881), 'k:')
hold off
xlim([0.5 9.3])
subplot(1,2,2)
bar(neurons_cc');
legend("GLM", "2nd Order GLM", "Staged Point-Process Model", "Position",[0.5  0.95  0  0], "Box","off", "Orientation","horizontal")
% savefig(h, 'results/final/models-DBR-cc-mPFC-8.fig')
% save("results/DBRccStatistical-mPFC-8.mat", "neurons_DBR", "neurons_cc")