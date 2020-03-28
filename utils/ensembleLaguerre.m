function X = ensembleLaguerre(x, H, lagnum, order, alpha)
%% laguerre basis
LagEx = zeros(H,lagnum);
for idx = 0:H-1
  for jdx =1:lagnum
    if idx<jdx
      tmps = sum((-1).^(0:idx).*arrayfun(@(x)nchoosek(idx,x), 0:idx).*arrayfun(@(x)nchoosek(jdx,x), 0:idx).*(1-alpha).^(0:idx).*alpha.^(idx:-1:0));
      LagEx(idx+1,jdx) = (-1)^idx*alpha^((jdx-idx)/2)*(1-alpha)^(1/2)*tmps;
    else
      tmps = sum((-1).^(0:jdx).*arrayfun(@(x)nchoosek(idx,x), 0:jdx).*arrayfun(@(x)nchoosek(jdx,x), 0:jdx).*(1-alpha).^(0:jdx).*alpha.^(jdx:-1:0));
      LagEx(idx+1,jdx) = (-1)^jdx*alpha^((idx-jdx)/2)*(1-alpha)^(1/2)*tmps;
    end
  end
end
%% For each channel, change the first order ensemble history to lag array
[~, channels] = size(x);
x_ensemble = ensemble(x, H);
[timeBins, ~] = size(x_ensemble);
x_new = zeros(timeBins, channels, lagnum);
% Train_new = zeros(validate_base - train_base, size(PMd, 1) ,lagnum);
% Validate_new = zeros(test_base - validate_base,size(PMd, 1),lagnum);
% Test_new = zeros(size(M1, 2)-opt.relevantTimeslot-test_base+1,size(PMd, 1),opt.lagnum);
for idx = 1:channels
  x_tmp = x_ensemble(:,idx:channels:channels*(H-1)+idx);
  x_new(:,idx,:) = x_tmp*LagEx;
end
%% construct specific order ensemble
Xzero  = ones(timeBins, 1); % zeros order
Xfirst = ones(timeBins, channels*lagnum); % first order
for idx=1:channels
  Xfirst(:,(idx-1)*lagnum+1:idx*lagnum) = squeeze(x_new(:,idx,:));
end

if (order==1)
  X = [Xfirst Xzero];
  return
end

Xsecond = ones(timeBins, (channels*lagnum+1)*channels*lagnum/2); % second order
unitTril = logical(tril(ones(channels*lagnum)));
for k=1:timeBins
  Xcross = x_new(k,:)' * x_new(k,:);
  Xsecond(k,:) = Xcross(unitTril)';
end

if (order==2)
  X = [Xsecond Xfirst Xzero];
  return
end
end