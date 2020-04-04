function lagX = ensembleLaguerreRecur(x, H, lagnum, order, alpha)
if (order~=2)
  disp("Not supported!")
  return;
end

[timeBins, channels] = size(x);
lagX = zeros(timeBins-H+1, channels*lagnum);

% j=0, tau=0
lagX(1,1:channels) = sqrt(1-alpha)*x(1,:);
% j>0, tau=0
if (lagnum>1)
  for idj=2:lagnum
    lagX(1,(idj-1)*channels+1:idj*channels) = sqrt(alpha^idj*(1-alpha))*x(1,:);
  end
end
% j=0, tau>0
for idx=2:timeBins-H+1
  lagX(idx, 1:channels) = sqrt(alpha)*lagX(idx-1,1:channels) + sqrt(1-alpha)*x(idx,:);
end
% j>0, tau>0
if (lagnum>1)
  for idx=2:2:timeBins-H+1
    for idj=2:lagnum
      lagX(idx, (idj-1)*channels+1:idj*channels) = ...
        sqrt(alpha)*lagX(idx-1, (idj-1)*channels+1:idj*channels) + ...
        sqrt(alpha)*lagX(idx, (idj-2)*channels+1:(idj-1)*channels) - ...
        lagX(idx-1, (idj-2)*channels+1:(idj-1)*channels);
    end
  end
end


end