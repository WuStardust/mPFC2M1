close all;clear;clc;
rng('default')
addpath models\
addpath utils\
load("data/data_rat010_0615_spike_train_selected.mat")

%% GLM
GLM_explore_H = struct( ...
  "H",{}, "xi",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, ...
  "W",{}, "L",{}, "DBR",{}, "Lval",{}, "LHistory",{} ...
  );
for H=1:10 % Temporal history. Explore within [1:10]
  xi=0.1; % Weights random initial range parameter
  alpha=0; % Regulization parameter
  thres = 1e-3; % stop error tolerance
  iterThres = 7; % stop after error over threshold $ times
  maxIter = 40; % max iteration num, over needs re-initial
  M1Idx = 1; % select M1 neuron
  splitFun = @(history)splitData(mPFCspike,M1spike(:,M1Idx),history); % choose splitData function as first Order
  verbose = 2;
  [W, L, DBR, Lval, LHistory] = runGLM(H, xi, thres, iterThres, maxIter, alpha, splitFun, verbose);
  GLM_explore_H(H) = struct( ...
    "H",H, "xi",xi, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, ...
    "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LHistory ...
    );
end
save("results\GLM_explore_H.mat", "GLM_explore_H")
%% 2nd Order GLM
GLM_sec_explore_H = struct( ...
  "H",{}, "xi",{}, "thres",{}, "iterThres",{}, ...
  "maxIter",{}, "alpha",{}, "M1Idx",{}, ...
  "W",{}, "L",{}, "DBR",{} ...
  );
for H=1:10
  % H = 1; % temporal history
  xi=0.1; % Weights random initial range parameter
  alpha=0; % Regulization parameter
  thres = 1e-3; % stop error tolerance
  iterThres = 7; % stop after error over threshold $ times
  maxIter = 40; % max iteration num, over needs re-initial
  M1Idx = 1; % select M1 neuron
  splitFun = @(history)splitDataSecOrder(mPFCspike,M1spike(:,M1Idx),history); % choose splitData function as second Order
  verbose = 2;
  [W, L, GLM_sec_DBR] = runGLM(H, xi, thres, iterThres, maxIter, alpha, splitFun, verbose);
  GLM_sec_explore_H(H) = struct( ...
    "H",H, "xi",xi, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, ...
    "W",W, "L",L, "DBR",GLM_sec_DBR ...
    );
end
save("results\GLM_sec_explore_H.mat", "GLM_sec_explore_H")
%% show results
% load data; uncoment following lines when needed
% close all;clear;clc;
% addpath models/
% load data/data_rat010_0615_spike_train_selected.mat
% load results/GLM_explore_H.mat
% load results/GLM_sec_explore_H.mat

% get DBRs from differnt model result
GLM_sec_DBR = zeros(1,10);
GLM_DBR = zeros(1,10);
for i=1:10
  GLM_DBR(i)     = GLM_explore_H(i).DBR;
  GLM_sec_DBR(i) = GLM_sec_explore_H(i).DBR;
end
% plot DBR-history result
figure(1)
plot(10:10:100, GLM_DBR, 'b')
hold on
plot(10:10:100, GLM_sec_DBR, 'g')
hold off
legend("GLM", "2nd-Order GLM")
xlabel("Length of mPFC history(msec)")
xlim([10 100])
ylabel("7th neuron DBR")
title("DBR-history")

H=8;
% get GLM params & test result
W     = GLM_explore_H(H).W;
M1Idx = GLM_explore_H(H).M1Idx;
[~,~,~,~,~,testX,~,~,testY] = splitData(mPFCspike,M1spike(:,M1Idx),H);
GLMtestLambdaYpre = GLMmodel(testX, W);
GLMccList = zeros(1, 1000/10);
% smooth result
for kernelSize=10:10:1000
  smoothedLambda    = gaussianSmooth(testY, kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMtestLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLMccList(kernelSize/10) = cc(2);
end
% 2nd Order GLM
W     = GLM_sec_explore_H(H).W;
M1Idx = GLM_sec_explore_H(H).M1Idx;
[~,~,~,~,~,testX,~,~,testY] = splitDataSecOrder(mPFCspike,M1spike(:,M1Idx),H);
GLMsectestLambdaYpre = GLMmodel(testX, W);
GLM_sec_ccList = zeros(1, 1000/10);
% smooth
for kernelSize=10:10:1000
  smoothedLambda    = gaussianSmooth(testY, kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMsectestLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLM_sec_ccList(kernelSize/10) = cc(2);
end

% optimal smooth kernel size result
kernelSize = 100;
smoothedLambda = gaussianSmooth(testY, kernelSize);
smoothedGLMLambdaPre = gaussianSmooth(GLMtestLambdaYpre, kernelSize);
smoothedGLMsecLambdaPre = gaussianSmooth(GLMsectestLambdaYpre, kernelSize);

% display results
spikeLength = length(testY);
index = 1:spikeLength;
t = index/100;
figure(2)
% real spike
subplot(4,1,2)
area(t, testY(index))
xlabel("time(sec)")
set(gca, 'TickLength', [0 0])
set(gca, 'ytick', [])
set(gca, 'box', 'off')
% smoothed results
subplot(4,1,3)
hold on
h{1} = plot(t, smoothedLambda(index), 'r');
h{2} = plot(t, smoothedGLMLambdaPre(index), 'b');
h{3} = plot(t, smoothedGLMsecLambdaPre(index), 'g');
hold off
xlabel("time(sec)")
ylabel("Firing rate")
% cc analyze
subplot(4,1,4)
hold on
plot(0.1:0.1:10, GLMccList, 'b')
plot(0.1:0.1:10, GLM_sec_ccList, 'g')
hold off
xlabel("kenel size(sec)")
xlim([0 10])
ylabel("CC")
title("CC-kernel size")
legend([h{1}; h{2}; h{3}], "Actual M1 spike", "GLM", "2nd-Order GLM", ...
  "Position",[0.5  0.95  0  0], "Box","off", "Orientation","horizontal")
subplot(4,1,1)
area(t, eventTrain(index))
figure(3)
hold on
h{1} = plot(t, smoothedLambda(index), 'r');
h{2} = plot(t, smoothedGLMLambdaPre(index), 'b');
h{3} = plot(t, smoothedGLMsecLambdaPre(index), 'g');
hold off
