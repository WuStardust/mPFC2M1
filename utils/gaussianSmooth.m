function smoothedSignal = gaussianSmooth(signal,kernelSize)
% Smooth the signal with gaussian kernel
% Get Gaussian kernel
kernel=normpdf(linspace(-3,3,kernelSize+1),0,1); % take values from 3 sigma region
kernel = kernel/sum(kernel);
smoothedSignal = conv(signal, kernel, 'same');
end
