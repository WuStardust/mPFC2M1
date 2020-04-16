function [xi, y, xAxis, B, DBR] = dbr(p, spikes)
% time scale transform
q = -log(1-p);
%% code by Shenghui
mask = spikes==1;
subs = cumsum(spikes);
subs(mask) = subs(mask) - 1;
subs = subs(find(spikes,1)+1:find(spikes,1,'last'));
qk    = q   (find(spikes,1)+1:find(spikes,1,'last'));
xi = accumarray(subs', qk', [], @(x) ...,
  sum(x(1:length(x)-1)) - log(...,
  1 - rand() * (...,
  1 - exp(-x(length(x)))...,
  )...,
  )...,
  );
%% code from https://github.com/iahncajigas/nSTAT/blob/master/Analysis.m
% spikeindicies=find(spikes==1);
% spikNum=length(spikeindicies);
% xi=zeros(spikNum-1,1);
% for r=1:spikNum-1
%   total = 0;
% 
%   ind1=spikeindicies(r);
%   ind2=spikeindicies(r+1);
% 
%   total=total+sum(q(ind1+1:ind2-1));
% 
%   delta=-(1/q(ind2))*log(1-rand()*(1-exp(-q(ind2))));
% 
%   total=total+q(ind2)*delta;
% 
%   xi(r)=total;
% end

y = 1 - exp(-xi);

% distance
B = sort(y);
spikNum=sum(spikes);
xAxis = 1/spikNum:1/spikNum:1-1/spikNum;
D = max(abs(B' - xAxis));
DBR = D / 1.36 * sqrt(spikNum);

%     plotDBR(y)
end
