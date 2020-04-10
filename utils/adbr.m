function avrDBR = adbr(p, spikes, spikeLength)
avrDBR = 0;
for i=1:100
  seg = 1;
  startIdx = 1;
  stopIdx  = 3000; % floor(spikeLength/5); % split to 10 segment for calculate
  step     = 1500; % floor(stopIdx/2); % 50 percent overlap
  DBR = 0;
  while(stopIdx<spikeLength)
    opt.DTCorrelation = 1;
    opt.sampleRate = 100;
    [~,~,~,~,temp]=computeKSStats(spikes(startIdx:stopIdx), p(startIdx:stopIdx)*100, opt);
    DBR     = DBR + temp/1.36*sqrt(sum(spikes(startIdx:stopIdx)));
    seg      = seg + 1;
    startIdx = step*(seg-1)+1;
    stopIdx  = step*(seg+1);
  end
  avrDBR = avrDBR + DBR/(seg-1);
end
avrDBR = avrDBR/100;
end
