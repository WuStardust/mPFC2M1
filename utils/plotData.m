function plotData(spikeTrainY, spikeTrainYpredict, lambdaYTrainPredict, LHistory, W)
    figure(1)

    t = 0:0.01:(length(spikeTrainY) - 1) * 0.01;
    
    subplot(5, 1, 1)
    plot(t, spikeTrainY);
    title('spikeY Ground Truth')

    subplot(5, 1, 2)
    plot(t, spikeTrainYpredict);
    title('spikeY predict')
    
    subplot(5, 1, 3);
    plot(t, lambdaYTrainPredict);
    hold on
    plot(t, gaussianSmooth(lambdaYTrainPredict,30), 'r')
    plot(t, gaussianSmooth(spikeTrainY,30), 'b')
    hold off
    legend("predict", "smoothedPridicted", "smoothedTruth")
    xlabel('Time(sec)')
    ylabel('lambdaY Predict')

    subplot(5, 1, 4)
    plot(LHistory)
    subplot(5, 1, 5)
    histogram(W)

    drawnow
end