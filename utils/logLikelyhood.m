function L = logLikelyhood(spikeTrainY, lambdaYTrainPredict, normW)
    L = spikeTrainY' * log(lambdaYTrainPredict) + (1 - spikeTrainY') * log(1 - lambdaYTrainPredict) - normW;

    if(isnan(L))
        disp('Error: L is NaN!');
    end
end
