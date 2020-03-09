function [W, L, DBR, Lval, LHistory] = runGLM(H, xi, threshold, iterationThres, maxIterations, alpha, splitFunc, verbose)
% Train & test GLM
%% Spike train ensemble & split train/test
[trainLen,~,testLen,trainX,valX,testX,trainY,valY,testY] = splitFunc(H);
%% Initialize params & variables
[~, Nx] = size(trainX);
W = xi / sqrt(Nx) * (2 * rand(1, Nx) - 1);

LHistory = zeros(1, maxIterations);
L = -Inf; % initialize Lpre as -Inf
overIterations = 0;
%% Train GLM
for iter=1:maxIterations
  % forward
  trainLambdaYpre = GLMmodel(trainX, W);
  % update Weigths
  G = (trainY - trainLambdaYpre)' * trainX - alpha * abs(W) ./ W;
  He = trainX' * spdiags(trainLambdaYpre.*(1-trainLambdaYpre),0,trainLen,trainLen) * trainX;
  if (rcond(He) < 1e-10)
    disp("Hessian matrix is singular. Re-initial or change hyper-param")
    break;
  end
  W = W + G / He;
  % validate
  [valLambdaYpre, valYpre] = GLMmodel(valX, W);
  Lnew = logLikelyhood(valY, valLambdaYpre, alpha*norm(W, 1));
  err = Lnew - L;
  if (err < threshold)
    overIterations = overIterations + 1;
  else
    overIterations = 0;
  end
  L = Lnew;
  LHistory(iter:length(LHistory)) = L; % record L
  if (overIterations > iterationThres)
    break;
  end
  if (verbose <= 0)
    disp([num2str(iter),'/',num2str(maxIterations),'...Lval=',num2str(L)]);
  end
end
if (verbose <= 1)
  plotData(valY(1:10000), valYpre(1:10000), valLambdaYpre(1:10000), LHistory, W)
end
if (verbose <= 2)
  disp([num2str(iter),'/',num2str(maxIterations),':Train Complete...L:',num2str(L), 9, '...H:', num2str(H), 9, '...xi:',num2str(xi), 9, '...alpha:',num2str(alpha)]);
end
Lval = L;
LHistory = LHistory(1:iter);
%% Test data
testLambdaYpre = GLMmodel(testX, W);
L = logLikelyhood(testY, testLambdaYpre, alpha*norm(W, 1));
DBR = adbr(testLambdaYpre, testY, testLen);
if (verbose <= 2)
  disp(['      Test Complete...L:',num2str(L), 9, '...DBR:',num2str(DBR)]);
end
end
