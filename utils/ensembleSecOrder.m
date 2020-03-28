function X = ensembleSecOrder(x, H)
%ensemble - get the signal second order ensemble with weight
%
% Syntax: X = ensemble(x, H)
%
% Get the signal ensemble from the input signal x and the history length H.
% For multiple channels input, the different channels will compact to signal channel
% Input: x -- [timeBins, channels] input signal
%        H -- history length
% Output: X -- [timeBins-H+1, channels*H+1] signal ensemble with weight
  [timeBins, channels] = size(x);
  X = zeros(timeBins-H+1, channels*H+1);
  for h=1:H
    X(:, channels*(h-1)+1:channels*h) = x(H-h+1:timeBins-h+1, :);
  end
  X(:, channels*H+1) = ones(timeBins-h+1, 1);
  
  XcovHat = zeros(timeBins-H+1, (channels*H-1)*channels*H/2);
  unitTril = logical(tril(ones(channels*H),-1)); % square of logical is same as before, so no need for square
  for k=H:timeBins
    Xhis = reshape(x(k-H+1:k,:)*1, channels*H, 1);
    Xcov = Xhis * Xhis';
    XcovHat(k-H+1,:) = Xcov(unitTril)';
  end

  X = [XcovHat, X];
end
