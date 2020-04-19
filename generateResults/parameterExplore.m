close all;clear;clc;
addpath models/
load data/data_rat010_0615_spike_train_selected_with_delay.mat
load results/ANN_explore_Nz.mat
load results/GLM_explore_H_1.mat
load results/GLM_sec_explore_H_1.mat
load results/ANN_explore_H_1.mat
ANN_explore_H2 = load("results\ANN_explore_H_1_Nz_6.mat", "ANN_explore_H");
ANN_explore_H2 = ANN_explore_H2.ANN_explore_H;

% Nz-DBR
ANN_Nz_DBR = zeros(1,25);
Nzlist = 1:25;
for Nz=Nzlist
  bestDBR = Inf;
  ANN_DBR_list = getFieldArray(ANN_explore_Nz, "DBR", Nz);
  ANN_Nz_DBR(Nz)=min(ANN_DBR_list);
end
h = figure("Name", "Nz-DBR");
plot(Nzlist, ANN_Nz_DBR(Nzlist))
xlabel("Number of hidden units")
ylabel("7th neuron DBR")
ylim([0.8,1])
title(['Nz-DBR(H=', num2str(ANN_explore_Nz(5,1).H), ')'])
savefig(h, "results/final/explore Nz-DBR.fig")

% H-DBR
% get DBRs from differnt model result
GLM_sec_DBR = zeros(1,20);
GLM_DBR = zeros(1,20);
ANN_DBR = zeros(1,20);
ANN2_DBR = zeros(1,20);
for H=2:20
  GLM_DBR(H)     = GLM_explore_H(H).DBR;
  GLM_sec_DBR(H) = GLM_sec_explore_H(H).DBR;
  ANN_DBR_list = getFieldArray(ANN_explore_H, "DBR", H);
  ANN_DBR(H) = min(ANN_DBR_list);
  ANN_DBR_list = getFieldArray(ANN_explore_H2, "DBR", H);
  ANN2_DBR(H) = min(ANN_DBR_list);
end
% plot DBR-history result
h = figure("Name", "history-DBR");
plot(20:10:200, GLM_DBR(2:20), 'b')
hold on
plot(20:10:200, GLM_sec_DBR(2:20), 'g')
plot(20:10:200, ANN_DBR(2:20), 'm')
plot(20:10:200, ANN2_DBR(2:20), '--m')
hold off
legend("GLM", "2nd-Order GLM", ['Staged Point-Process(Nz=', num2str(ANN_explore_H(2,1).Nz), ')'], ['Staged Point-Process(Nz=', num2str(ANN_explore_H2(2,1).Nz), ')'])
xlabel("Length of mPFC history(msec)")
xlim([10 200])
ylabel("7th neuron DBR")
title("history-DBR")
savefig(h, "results/final/explore history-DBR.fig")