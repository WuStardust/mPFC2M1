addpath utils/
load data/data_rat010_0615_spike_train_top_9.mat

t = 1565:0.01:1595;
index = 156500:159500;
h = figure("Name", "records");
subplot(4,1,1)
area(t, segTrain(index))
% xlabel("time(sec)")
ylabel("behavior")
subplot(4,1,2)
area(t, M1spike(index, 3), 'FaceAlpha',.3)
hold on
plot(t, gaussianSmooth(M1spike(index, 3), 200)*5, 'LineWidth',2)
hold off
% xlabel("time(sec)")
ylabel("7th M1")
set(gca, 'TickLength', [0 0])
set(gca, 'ytick', [])
set(gca, 'box', 'off')
subplot(4,1,3)
area(t, mPFCspike(index, 7), 'FaceAlpha',.3)
hold on
plot(t, gaussianSmooth(mPFCspike(index, 7), 200)*5, 'LineWidth',2)
hold off
% xlabel("time(sec)")
ylabel("7th mPFC")
set(gca, 'TickLength', [0 0])
set(gca, 'ytick', [])
set(gca, 'box', 'off')
subplot(4,1,4)
area(t, mPFCspike(index, 8), 'FaceAlpha',.3)
hold on
plot(t, gaussianSmooth(mPFCspike(index, 8), 200)*5, 'LineWidth',2)
hold off
xlabel("time(sec)")
ylabel("8th mPFC")
set(gca, 'TickLength', [0 0])
set(gca, 'ytick', [])
set(gca, 'box', 'off')

savefig(h, 'results/final/recordsEg.fig')
