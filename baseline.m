close all;clear;clc;
rng('default')
addpath models\
addpath utils\
load("data/data_rat010_0615_spike_train_selected_with_delay.mat")
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
  iterThres = 0; % stop after error over threshold $ times
  maxIter = 40; % max iteration num, over needs re-initial
  M1Idx = 1; % select M1 neuron
  splitFun = @(history)splitDataAdvance(1,mPFCspike,M1spike(:,M1Idx),eventTrain,optimalDelay(M1Idx),segTrain,history); % choose splitData function as first Order
  verbose = 2;
  [W, L, DBR, Lval, LHistory] = runGLM(H, xi, thres, iterThres, maxIter, alpha, splitFun, verbose);
  GLM_explore_H(H) = struct( ...
    "H",H, "xi",xi, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, ...
    "W",W, "L",L, "DBR",DBR, "Lval",Lval, "LHistory",LHistory ...
    );
end
save("results\GLM_explore_H_1.mat", "GLM_explore_H")
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
  iterThres = 0; % stop after error over threshold $ times
  maxIter = 40; % max iteration num, over needs re-initial
  M1Idx = 1; % select M1 neuron
  splitFun = @(history)splitDataAdvance(2,mPFCspike,M1spike(:,M1Idx),eventTrain,optimalDelay(M1Idx),segTrain,history);
  verbose = 2;
  [W, L, GLM_sec_DBR] = runGLM(H, xi, thres, iterThres, maxIter, alpha, splitFun, verbose);
  GLM_sec_explore_H(H) = struct( ...
    "H",H, "xi",xi, "thres",thres, "iterThres",iterThres, ...
    "maxIter",maxIter, "alpha",alpha, "M1Idx",M1Idx, ...
    "W",W, "L",L, "DBR",GLM_sec_DBR ...
    );
end
save("results\GLM_sec_explore_H_1.mat", "GLM_sec_explore_H")
%% show results
% load data; uncoment following lines when needed
% close all;clear;clc;
% addpath models/
% load data/data_rat010_0615_spike_train_selected_with_delay.mat
% load results/GLM_explore_H_new.mat
% load results/GLM_sec_explore_H_new.mat

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

H=10;
% get GLM params & test result
W     = GLM_explore_H(H).W;
M1Idx = GLM_explore_H(H).M1Idx;
[~,~,~,trainX,~,testX,trainY,~,testY] = splitDataAdvance(1,mPFCspike,M1spike(:,M1Idx),eventTrain,optimalDelay(M1Idx),segTrain,H);
GLMtestLambdaYpre = GLMmodel(testX, W);
GLMtrainLambdaYpre = GLMmodel(trainX,W);
GLMccListTest = zeros(1, 1000/10);
GLMccListTrain = zeros(1, 1000/10);
% smooth result
for kernelSize=10:10:1000
  smoothedLambda    = gaussianSmooth(testY, kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMtestLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLMccListTest(kernelSize/10) = cc(2);
  smoothedLambda    = gaussianSmooth(trainY, kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMtrainLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLMccListTrain(kernelSize/10) = cc(2);
end
% 2nd Order GLM
W     = GLM_sec_explore_H(H).W;
M1Idx = GLM_sec_explore_H(H).M1Idx;
[~,~,~,trainX,~,testX,trainY,~,testY,trainEvent,~,testEvent] = splitDataAdvance(2,mPFCspike,M1spike(:,M1Idx),eventTrain,optimalDelay(M1Idx),segTrain,H);
GLMsectestLambdaYpre = GLMmodel(testX, W);
GLMsectrainLambdaYpre = GLMmodel(trainX, W);
GLM_sec_ccListTest = zeros(1, 1000/10);
GLM_sec_ccListTrain = zeros(1, 1000/10);
% smooth
for kernelSize=10:10:1000
  smoothedLambda    = gaussianSmooth(testY, kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMsectestLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLM_sec_ccListTest(kernelSize/10) = cc(2);
  smoothedLambda    = gaussianSmooth(trainY, kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMsectrainLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLM_sec_ccListTrain(kernelSize/10) = cc(2);
end

% optimal smooth kernel size result
kernelSize = 30;
smoothedLambda = gaussianSmooth(trainY, kernelSize);
smoothedGLMLambdaPre = gaussianSmooth(GLMtrainLambdaYpre, kernelSize);
smoothedGLMsecLambdaPre = gaussianSmooth(GLMsectrainLambdaYpre, kernelSize);

% display results
% spikeLength = length(trainY);
% index = 1:spikeLength;
% t = index/100;
t = 55:0.01:100;
index = 5500:10000;
figure("Name", "Train")
subplot(4,1,1)
% area((1:length(trainEvent))/100, trainEvent)
area(t, trainEvent(index))
% real spike
subplot(4,1,2)
area(t, trainY(index))
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
plot(0.1:0.1:10, GLMccListTrain, 'b')
plot(0.1:0.1:10, GLM_sec_ccListTrain, 'g')
hold off
xlabel("kenel size(sec)")
xlim([0 10])
ylabel("CC")
title("CC-kernel size")
legend([h{1}; h{2}; h{3}], "Actual M1 spike", "GLM", "2nd-Order GLM", ...
  "Position",[0.5  0.95  0  0], "Box","off", "Orientation","horizontal")

% optimal smooth kernel size result
kernelSize = 30;
smoothedLambda = gaussianSmooth(testY, kernelSize);
smoothedGLMLambdaPre = gaussianSmooth(GLMtestLambdaYpre, kernelSize);
smoothedGLMsecLambdaPre = gaussianSmooth(GLMsectestLambdaYpre, kernelSize);

% display results
% spikeLength = length(testY);
% index = 1:spikeLength;
% t = index/100;
t = 200:0.01:300;
index = 20000:30000;
figure("Name", "Test")
subplot(4,1,1)
% area((1:length(testEvent))/100, testEvent)
area(t, testEvent(index))
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
plot(0.1:0.1:10, GLMccListTest, 'b')
plot(0.1:0.1:10, GLM_sec_ccListTest, 'g')
hold off
xlabel("kenel size(sec)")
xlim([0 10])
ylabel("CC")
title("CC-kernel size")
legend([h{1}; h{2}; h{3}], "Actual M1 spike", "GLM", "2nd-Order GLM", ...
  "Position",[0.5  0.95  0  0], "Box","off", "Orientation","horizontal")
