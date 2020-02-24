function avrDBR = adbr(p, spikes, spikeLength)
avrDBR = 0;
for i=1:100
  seg = 1;
  startIdx = 1;
  stopIdx  = 3000; % floor(spikeLength/5); % split to 10 segment for calculate
  step     = 1500; % floor(stopIdx/2); % 50 percent overlap
  DBR = 0;
  while(stopIdx<spikeLength)
    DBR      = DBR + dbr(p(startIdx:stopIdx)', spikes(startIdx:stopIdx)');
    seg      = seg + 1;
    startIdx = step*(seg-1)+1;
    stopIdx  = step*(seg+1);
  end
  avrDBR = avrDBR + DBR/(seg-1);
end
avrDBR = avrDBR/100;
end
