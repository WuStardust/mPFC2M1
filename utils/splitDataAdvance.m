function [trainLen, valLen, testLen, trainX, valX, testX, trainY, valY, testY, trainEvent, valEvent, testEvent] = splitDataAdvance(order, inputSignal, outputSignal, eventTrain, delay, segTrain, H)
%% get data ensemble with history
pattern = 'vm(\d+)-(\d+)-(\d+[.]\d+)';
if (order==1)
  Xhat = ensemble(inputSignal, H);
elseif (order==2)
  Xhat = ensembleSecOrder(inputSignal, H);
elseif (regexp(order,pattern))
  regOut = regexp(order,pattern,'tokens');
  lagOrder = str2double(regOut{1}{1});
  lagnum   = str2double(regOut{1}{2});
  lagalpha = str2double(regOut{1}{3});
  Xhat = ensembleLaguerre(inputSignal, H, lagnum, lagOrder, lagalpha);
end

%% success trail number
successCount = 0;
for i=2:length(eventTrain)
  if (eventTrain(i-1) == 1 && eventTrain(i) == 0)
    successCount = successCount + 1;
  end
end
%% train set
trailNo = 1;
trainLen = 0;
trainX = zeros(size(Xhat));
trainY = zeros(size(outputSignal));
trainEvent = zeros(size(eventTrain));
trainNum = ceil(successCount*0.6);
for i=H+delay:length(eventTrain)
  if (eventTrain(i-1)==0 && eventTrain(i) == 1)
    trailNo = trailNo + 1;
  end
  
  isTrainTrail = (trailNo <= trainNum);
  if (isTrainTrail && segTrain(i) > 0)
    trainLen = trainLen + 1;
    trainX(trainLen,:) = Xhat(i-H-delay+1,:);
    trainY(trainLen) = outputSignal(i-delay);
    trainEvent(trainLen) = segTrain(i);
  end
  
  if (~isTrainTrail); break; end
end
trainX = trainX(1:trainLen,:);
trainY = trainY(1:trainLen);
trainEvent = trainEvent(1:trainLen);
%% validate & test
valStart = i-H-delay+1;
valLen = floor((length(Xhat)-valStart)/2);
valX = Xhat(valStart+1:valStart+valLen,:);
valY = outputSignal(valStart+H:valStart+valLen+H-1);
valEvent = segTrain(valStart+H+delay:valStart+valLen+H+delay-1);

testX = Xhat(valStart+valLen+1:end,:);
testY = outputSignal(valStart+valLen+H:end);
testLen = length(testY);
testEvent = segTrain(valStart+valLen+H+delay:end);
end