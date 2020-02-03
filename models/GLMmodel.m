function [lambdaYpredict, spikeYpredict] = GLMmodel(spikeTrainEnsemble, W)
% GLM model
precise = 1e-6; % deal with tructed error of matlab
lambdaYpredict = sigmaFunc(spikeTrainEnsemble*W');
lambdaYpredict(lambdaYpredict<precise) = precise;
lambdaYpredict(1-lambdaYpredict<precise) = 1-precise;
spikeYpredict = lambda2Spike(lambdaYpredict);
end