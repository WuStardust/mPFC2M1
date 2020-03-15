function [M1, mPFC, ratios, eventTrain, M1num, mPFCnum, actTrain, segTrain, segRatios] = spikeTime2Train(filename)
S = load(filename);
timebins  = 0.01; % 10ms time bins
spikeLength = 1e7;
M1num       = 0;
mPFCnum     = 0;
%% get all channel names
for channelNo = 1:32
  for sub = ['a', 'b']
    channelName = ['WB', num2str(channelNo, '%02d'), sub];
    if (~isfield(S, channelName))
      continue;
    end
    if (channelNo <= 16)
      M1num     = M1num + 1;
    else
      mPFCnum   = mPFCnum + 1;
    end
    spikeLength = min(spikeLength, ceil(max(eval(['S.', channelName]))/timebins));
  end
end
M1   = zeros(spikeLength, M1num);
mPFC = zeros(spikeLength, mPFCnum);

%% from spike time to spike train
M1i = 1; mPFCi = 1;
for channelNo = 1:32
  for sub = ['a', 'b']
    channelName = ['WB', num2str(channelNo, '%02d'), sub];
    if (~isfield(S, channelName))
      continue;
    end
    
    spikeTrain = accumarray(fix(eval(['S.', channelName])*1/timebins)+1, 1);
    if (channelNo <= 16) % first 16 is M1
      M1(:, M1i) = spikeTrain(1:spikeLength);
      M1i = M1i+1;
    else % latter 16 is mPFC
      mPFC(:,mPFCi) = spikeTrain(1:spikeLength);
      mPFCi = mPFCi+1;
    end
  end
end
%% event time to event train
% highLever = accumarray(fix(S.EVT01(:,1))/timebins+1, 1);
% lowLever  = accumarray(fix(S.EVT02(:,1))/timebins+1, 1);
success   = accumarray(fix(S.EVT03(:,1)/timebins)+1, 1)>0;
startTrail= accumarray(fix(S.EVT05(:,1)/timebins)+1, 1)>0;
press     = accumarray(fix(S.EVT06(:,1))/timebins+1, 1);
release   = accumarray(fix(S.EVT07(:,1)/timebins)+1, 1)>0;

suc = false;
eventTrain = zeros(1, spikeLength);
segTrain = zeros(1,spikeLength);
n=0; p=0;
for i=1:spikeLength
  if (i>length(startTrail)); continue; end
  if (startTrail(i) == 1)
    startTime = i;
    suc = false;
    continue;
  end
  
  if (i>length(press)); continue; end
  if (press(i) == 1)
    pressTime = i;
    continue;
  end
  
  if (i>length(success)); continue; end
  if (success(i) == 1)
    suc = true;
    p=p+1;
  end
  
  if (i>length(release)); continue; end
  if (release(i) == 1 && suc == true)
    eventTrain(startTime:i) = 1;
    segTrain(startTime:startTime+50) = 1;
    segTrain(pressTime:pressTime+50) = 2;
    segTrain(i:i+50) = 3;
    n=n+1;
    suc = false;
  end
end

actTrain = zeros(1,spikeLength);
for i=1:spikeLength
  if (i>length(press)); continue; end
  if (press(i) == 1)
    presstime = i;
    actTrain(max(i-10,1):i+10) = 1;
  end

  if (i>length(release)); continue; end
  if (release(i) == 1)
    actTrain(presstime:i) = 1;
    actTrain(max(i-10,1):i+10) = 1;
  end
end
%% calculate spike/total ratio & multi-spike/spike ratio
% all data
M1withSpike = mean(sum(M1>0)/spikeLength);
mPFCwithSpike = mean(sum(mPFC>0)/spikeLength);
M1multiSpike = mean(sum(M1>1)./sum(M1>0));
mPFCmultiSpike = mean(sum(mPFC>1)./sum(mPFC>0));
ratios = [M1withSpike, mPFCwithSpike, M1multiSpike, mPFCmultiSpike];
% seg data
segM1withSpike = mean(sum(M1.*(segTrain'>0)>0)/sum(segTrain>0));
segmPFCwithSpike = mean(sum(mPFC.*(segTrain'>0)>0)/sum(segTrain>0));
segM1multiSpike = mean(sum(M1.*(segTrain'>0)>1)./sum(M1.*(segTrain'>0)>0));
segmPFCmultiSpike = mean(sum(mPFC.*(segTrain'>0)>1)./sum(mPFC.*(segTrain'>0)>0));
segRatios = [segM1withSpike, segmPFCwithSpike, segM1multiSpike, segmPFCmultiSpike];
% mutispike also represented by 1
M1 = double(M1>0);
mPFC = double(mPFC>0);
end
