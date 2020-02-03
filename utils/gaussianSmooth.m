function smoothedSignal = gaussianSmooth(signal,kernelSize)
% Smooth the signal with gaussian kernel
%   The Gaussian kernel has 0 mean, unit stander deviation, specified
%   kernel size.
%   ATTENTION: unit is second & 10msec per point
% Get Gaussian kernel
kernel = exp(-((-kernelSize:kernelSize)/100).^2) / sqrt(2*pi);
kernel = kernel/sum(kernel);
smoothedSignal = conv(signal, kernel, 'same');
end
