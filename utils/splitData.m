function [trainLen, valLen, testLen, trainX, valX, testX, trainY, valY, testY] = splitData(inputSignal,outputSignal,H)
%Split data into train, val, test
Xhat = ensemble(inputSignal, H);
trainLen = ceil(0.6*length(inputSignal));
valLen   = ceil(0.2*length(inputSignal));
testLen  = length(Xhat) - trainLen - valLen;

trainX = Xhat(1:trainLen,:);
valX   = Xhat(trainLen+1:trainLen+valLen,:);
testX  = Xhat(trainLen+valLen+1:end,:);
trainY = outputSignal(H:H+trainLen-1);
valY   = outputSignal(H+trainLen:H+trainLen+valLen-1);
testY  = outputSignal(H+trainLen+valLen:end);
end
