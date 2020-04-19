close all;clear;clc;
addpath models/
load data/data_rat010_0615_spike_train_selected_with_delay.mat
load results/GLM_explore_H_1.mat
load results/GLM_sec_explore_H_1.mat
load results/ANN_explore_H_1_Nz_6.mat
%% prediction
H=20;
trainT = 55:0.01:85;
trainIndex = 5500:8500;
testT = 170:0.01:200;
testIndex = 17000:20000;
% get GLM params & test result
W     = GLM_explore_H(H).W;
M1Idx = GLM_explore_H(H).M1Idx;
[~,~,~,trainX,~,testX,trainY,~,testY] = splitDataAdvance(1,mPFCspike,M1spike(:,M1Idx),eventTrain,optimalDelay(M1Idx),segTrain,H);
GLMtestLambdaYpre = GLMmodel(testX(testIndex,:), W);
GLMtrainLambdaYpre = GLMmodel(trainX(trainIndex,:),W);
GLMccListTest = zeros(1, 1000/10);
GLMccListTrain = zeros(1, 1000/10);
% smooth result
for kernelSize=10:10:1000
  smoothedLambda    = gaussianSmooth(testY(testIndex), kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMtestLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLMccListTest(kernelSize/10) = cc(2);
  smoothedLambda    = gaussianSmooth(trainY(trainIndex), kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMtrainLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLMccListTrain(kernelSize/10) = cc(2);
end
% 2nd Order GLM
W     = GLM_sec_explore_H(H).W;
M1Idx = GLM_sec_explore_H(H).M1Idx;
lagnum = GLM_sec_explore_H(H).lagnum;
lagalpha = GLM_sec_explore_H(H).lagalpha;
order = ['vm2-',num2str(lagnum),'-',num2str(lagalpha)];
[~,~,~,trainX,~,testX,trainY,~,testY,trainEvent,~,testEvent] = splitDataAdvance(order,mPFCspike,M1spike(:,M1Idx),eventTrain,optimalDelay(M1Idx),segTrain,H);
GLMsectestLambdaYpre = GLMmodel(testX(testIndex,:), W);
GLMsectrainLambdaYpre = GLMmodel(trainX(trainIndex,:), W);
GLM_sec_ccListTest = zeros(1, 1000/10);
GLM_sec_ccListTrain = zeros(1, 1000/10);
% smooth
for kernelSize=10:10:1000
  smoothedLambda    = gaussianSmooth(testY(testIndex), kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMsectestLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLM_sec_ccListTest(kernelSize/10) = cc(2);
  smoothedLambda    = gaussianSmooth(trainY(trainIndex), kernelSize);
  smoothedLambdaPre = gaussianSmooth(GLMsectrainLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  GLM_sec_ccListTrain(kernelSize/10) = cc(2);
end
% ANN
tmp = getFieldArray(ANN_explore_H, "DBR", H);
[~,idx] = min(tmp);
W     = ANN_explore_H(H,idx).W;
M1Idx = ANN_explore_H(H,idx).M1Idx;
Nz    = ANN_explore_H(H,idx).Nz;
[~,~,~,trainX,~,testX,trainY,~,testY] = splitDataAdvance(1,mPFCspike,M1spike(:,M1Idx),eventTrain,optimalDelay(M1Idx),segTrain,H);
[~, Nx] = size(testX);
ANNtestLambdaYpre  = ANNmodel(testX(testIndex,:),  W, Nx, Nz);
ANNtrainLambdaYpre = ANNmodel(trainX(trainIndex,:), W, Nx, Nz);
ANNccListTest = zeros(1, 1000/10);
ANNccListTrain = zeros(1, 1000/10);
% smooth
for kernelSize=10:10:1000
  smoothedLambda    = gaussianSmooth(testY(testIndex), kernelSize);
  smoothedLambdaPre = gaussianSmooth(ANNtestLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  ANNccListTest(kernelSize/10) = cc(2);
  smoothedLambda    = gaussianSmooth(trainY(trainIndex), kernelSize);
  smoothedLambdaPre = gaussianSmooth(ANNtrainLambdaYpre, kernelSize);
  cc = corrcoef(smoothedLambda, smoothedLambdaPre);
  ANNccListTrain(kernelSize/10) = cc(2);
end
%% train
% optimal smooth kernel size result
kernelSize = 100;
smoothedLambda = gaussianSmooth(trainY(trainIndex), kernelSize);
smoothedGLMLambdaPre = gaussianSmooth(GLMtrainLambdaYpre, kernelSize);
smoothedGLMsecLambdaPre = gaussianSmooth(GLMsectrainLambdaYpre, kernelSize);
smoothedANNLambdaPre = gaussianSmooth(ANNtrainLambdaYpre, kernelSize);

% display results
% spikeLength = length(trainY);
% index = 1:spikeLength;
% t = index/100;
% t = 55:0.01:85;
% index = 5500:8500;
h = figure("Name", "Train");
subplot(4,1,1)
% area((1:length(trainEvent))/100, trainEvent)
area(trainT, trainEvent(trainIndex))
% real spike
subplot(4,1,2)
area(trainT, trainY(trainIndex))
xlabel("time(sec)")
set(gca, 'TickLength', [0 0])
set(gca, 'ytick', [])
set(gca, 'box', 'off')
% smoothed results
subplot(4,1,3)
l{1} = plot(trainT, smoothedLambda, 'r');
hold on
l{2} = plot(trainT, smoothedGLMLambdaPre, 'b');
l{3} = plot(trainT, smoothedGLMsecLambdaPre, 'g');
l{4} = plot(trainT, smoothedANNLambdaPre, 'm');
hold off
xlabel("time(sec)")
ylabel("Firing rate")
% cc analyze
subplot(4,1,4)
plot(0.1:0.1:10, GLMccListTrain, 'b')
hold on
plot(0.1:0.1:10, GLM_sec_ccListTrain, 'g')
plot(0.1:0.1:10, ANNccListTrain, 'm');
hold off
xlabel("kenel size(sec)")
xlim([0 10])
ylabel("CC")
title("CC-kernel size")
legend([l{1}; l{2}; l{3}; l{4}], "Actual M1 spike", "GLM", "2nd-Order GLM", ...
  "Staged Point-Process Mode", "Position",[0.5  0.95  0  0], ...
  "Box","off", "Orientation","horizontal")
savefig(h, ['results/final/trainEg-neuron-7-Nz-', num2str(Nz), '-H-', num2str(H), '.fig'])
%% test
% optimal smooth kernel size result
kernelSize = 100;
smoothedLambda = gaussianSmooth(testY(testIndex), kernelSize);
smoothedGLMLambdaPre = gaussianSmooth(GLMtestLambdaYpre, kernelSize);
smoothedGLMsecLambdaPre = gaussianSmooth(GLMsectestLambdaYpre, kernelSize);
smoothedANNLambdaPre = gaussianSmooth(ANNtestLambdaYpre, kernelSize);

% display results
% spikeLength = length(testY);
% index = 1:spikeLength;
% t = index/100;
% t = 170:0.01:200;
% index = 17000:20000;
h = figure("Name", "Test");
subplot(4,1,1)
% area((1:length(testEvent))/100, testEvent)
area(testT, testEvent(testIndex))
% real spike
subplot(4,1,2)
area(testT, testY(testIndex))
xlabel("time(sec)")
set(gca, 'TickLength', [0 0])
set(gca, 'ytick', [])
set(gca, 'box', 'off')
% smoothed results
subplot(4,1,3)
l{1} = plot(testT, smoothedLambda, 'r');
hold on
l{2} = plot(testT, smoothedGLMLambdaPre, 'b');
l{3} = plot(testT, smoothedGLMsecLambdaPre, 'g');
l{4} = plot(testT, smoothedANNLambdaPre, 'm');
xlabel("time(sec)")
ylabel("Firing rate")
% cc analyze
subplot(4,1,4)
plot(0.1:0.1:10, GLMccListTest, 'b')
hold on
plot(0.1:0.1:10, GLM_sec_ccListTest, 'g')
plot(0.1:0.1:10, ANNccListTest, 'm');
hold off
xlabel("kenel size(sec)")
xlim([0 10])
ylabel("CC")
title("CC-kernel size")
legend([l{1}; l{2}; l{3}; l{4}], "Actual M1 spike", "GLM", "2nd-Order GLM", ...
  "Staged Point-Process Model", "Position",[0.5  0.95  0  0], ...
  "Box","off", "Orientation","horizontal")
savefig(h, ['results/final/testEg-neuron-7-Nz-', num2str(Nz), '-H-', num2str(H), '.fig'])
%% DBR
opt.DTCorrelation = 1;
opt.sampleRate = 100;
h = figure("Name","ks plot");
hold on
[~,~, xAxis, KSSorted]=computeKSStats(testY(testIndex), GLMtestLambdaYpre*100, opt);
plot(xAxis, KSSorted, 'b')
[~,~, xAxis, KSSorted]=computeKSStats(testY(testIndex), GLMsectestLambdaYpre*100, opt);
plot(xAxis, KSSorted, 'g')
[~,~, xAxis, KSSorted]=computeKSStats(testY(testIndex), ANNtestLambdaYpre*100, opt);
plot(xAxis, KSSorted, 'm')
plot(xAxis, xAxis, 'k')
bound = 1.36/sqrt(sum(testY(testIndex)));
plot(bound:0.001:1, 0:0.001:1-bound, 'k:')
plot(0:0.001:1-bound, bound:0.001:1, 'k:')
hold off
legend("GLM", "2nd-Order GLM", "Staged Point-Process Model", "Position",[0.5  0.95  0  0], "Box","off", "Orientation","horizontal")
% title("neuron 7 ks plot")
savefig(h, "results/final/DBREg.fig")