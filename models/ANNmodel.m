function [lambdaYpredict, spikeYpredict, lambdaZ] = ANNmodel(spikeTrainEnsemble, W, Nx, Nz)
% ANN model
precise = 1e-6; % deal with tructed error of matlab
w = reshape(W(1:Nx*Nz), Nx, Nz)';
theta = W(Nx*Nz+1:length(W));

% hidden layer
lambdaZ = sigmaFunc(spikeTrainEnsemble*w');
lambdaZ(lambdaZ<precise) = precise;
lambdaZ(1-lambdaZ<precise) = 1-precise;

% output layer
lambdaYpredict = sigmaFunc([lambdaZ, ones(length(lambdaZ), 1)]*theta');
lambdaYpredict(lambdaYpredict<precise) = precise;
lambdaYpredict(1-lambdaYpredict<precise) = 1-precise;

spikeYpredict = lambda2Spike(lambdaYpredict);
end
