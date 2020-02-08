function [W,L,DBR,Whis] = runANN(H, Nz, xi1, xi2, mu, threshold, ...
  iterationThres, maxIterations, alpha, splitFunc, verbose)
% Train & test ANN
%% Spike train ensemble & split train/test
[trainLen,valLen,testLen,trainX,valX,testX,trainY,valY,testY] = splitFunc(H);
%% Initialize params & variables
[~, Nx] = size(trainX);
w     = xi2 / sqrt(Nx) * (2*rand(Nx, Nz  ) - 1);
theta = xi1 / sqrt(Nz) * (2*rand(1 , Nz+1) - 1);
W = [reshape(w, 1, Nx*Nz), theta];

LtrainHis = zeros(1, maxIterations);
LvalHis = zeros(1, maxIterations);
iter = 1;
valConverge = 0;
Whis = zeros(maxIterations+1,length(W));
Whis(1,:) = W;

%% train ANN
% forward
[trainLambdaYpre, ~, trainLambdaZpre] = ANNmodel(trainX, W, H, Nx, Nz);
Ltrain = logLikelyhood(trainY, trainLambdaYpre, alpha*norm(W, 1));
LtrainHis(iter:end) = Ltrain;
valLambdaYpre= ANNmodel(valX, W, H, Nx, Nz);
Lval = logLikelyhood(valY, valLambdaYpre, alpha*norm(W, 1));
LvalHis(iter:end) = Lval;
LvalBest = Lval;
bestIter = iter;

seg = 1;
startIdx = 1;
stopIdx  = 10000;
DBRval = 0;
while(stopIdx<valLen)
  DBRval   = DBRval + dbr(valLambdaYpre(startIdx:stopIdx)', valY(startIdx:stopIdx)');
  seg      = seg + 1;
  startIdx = 5000*(seg-1)+1;
  stopIdx  = 5000*(seg+1);
end
DBRval = DBRval/(seg-1);

while (iter<maxIterations)
  if (verbose<=0)
    disp(strcat(num2str(iter-1),'/',num2str(maxIterations), ...
      '...L=',num2str(Ltrain),'...mu=',num2str(mu),'...Lval=',num2str(Lval), ...
      '...DBRval',num2str(DBRval)));
  end

  % BP
  theta = W(Nx*Nz+1:Nx*Nz+Nz);
  % ----- Gradient
  dSpikeLambdaY = trainY - trainLambdaYpre;
  dSpikeLambdaYZProduct = dSpikeLambdaY .* trainLambdaZpre.* (1 - trainLambdaZpre);
  thetaDProdcut = theta .* dSpikeLambdaYZProduct;

  Gtheta = dSpikeLambdaY' * [trainLambdaZpre, ones(trainLen, 1)];
  
  Gw = zeros(1, Nx*Nz);
  for h=1:Nz
    Gw((h-1)*Nx+1:h*Nx) = thetaDProdcut(:, h)'*trainX;
  end

  Grad = [Gw, Gtheta] - alpha * abs(W)./W;
  % ----- Hessian
  LambdaYPro = (1 - trainLambdaYpre) .* trainLambdaYpre;
  
  Hetheta2 = - [trainLambdaZpre, ones(trainLen, 1)]' * ...
    spdiags(LambdaYPro,0,trainLen,trainLen) * [trainLambdaZpre, ones(trainLen, 1)];
  
  Hew2 = zeros(Nz * Nx);
  Hewtheta = zeros(Nz+1, Nz*Nx);
  
  LambdaZPro = trainLambdaZpre .* (1 - trainLambdaZpre);
  thetaLambdaZPro = theta .* LambdaZPro;
  thetaYZPro = LambdaYPro .* thetaLambdaZPro;
  
  bias = thetaLambdaZPro .* dSpikeLambdaY .* (1 - 2 * trainLambdaZpre);
  wthetaBias = dSpikeLambdaY .* LambdaZPro;
  biasM = eye(Nz);
  for m=1:Nz
    diagV = - thetaYZPro(:,m).*thetaLambdaZPro + biasM(m,:).*bias;
    for n=1:Nz
      Hew2((m-1)*Nx+1:m*Nx, (n-1)*Nx+1:n*Nx) = trainX' * ...
        spdiags(diagV(:, n),0,trainLen,trainLen) * trainX;
    end
    Hewtheta(1:Nz, Nx*(m-1)+1:Nx*m) = (- thetaYZPro(:,m).*trainLambdaZpre + ...
      biasM(m,:).*wthetaBias)' * trainX;
  end
  Hewtheta(Nz+1, :) = reshape((- thetaYZPro' * trainX)', 1, Nz*Nx);
  
  Hessian = [Hew2, Hewtheta'; Hewtheta, Hetheta2];
  if (rcond(Hessian) < 1e-15 || sum(sum(isnan(Hessian))) == 1)
    disp('Error: BAD He!');
  end
  
  while(1)
    % lm update
    Wnew = W - Grad/(Hessian-mu*eye(size(Hessian)));

    % forward
    [trainLambdaYpre, ~, trainLambdaZpre] = ANNmodel(trainX, Wnew, H, Nx, Nz);
    Lnew = logLikelyhood(trainY, trainLambdaYpre, alpha*norm(Wnew, 1));
    
    % check L
    if (Lnew<Ltrain) % fail to get better L, use lm with better mu.
      if (mu < 1e7)
        % update mu
        mu = mu*1000;
        continue; % increase mu and do it again
      else
        break; % mu is too large, still go on
      end
    else
      break; % Lnew is better, nice to go on
    end
  end
  
  % step: Lnew is better OR mu is too Large
  iter = iter + 1;
  W = Wnew;
  Whis(iter,:) = W;
  LtrainHis(iter) = Lnew;
  % update mu
  if (Lnew>Ltrain) % for Lnew is smaller but mu is large, don't change mu
    if (mu>1e-7) % mu is not so small, change mu; if mu is small, do not change
      mu = mu/100;
    end
  end
  Ltrain = Lnew;
  
  % validate
  [valLambdaYpre, valYpre] = ANNmodel(valX, W, H, Nx, Nz);
  LvalNew = logLikelyhood(valY, valLambdaYpre, alpha*norm(W, 1));
  if (LvalBest <= LvalNew)
    LvalBest = LvalNew;
    bestIter = iter;
  end
  LvalBest = max(LvalBest, LvalNew);
  % L on validation set change too little, or drop too much
  if (abs(LvalNew-Lval)<threshold || LvalNew-LvalBest<-50)
    valConverge = valConverge + 1;
  else
    valConverge = 0;
  end
  Lval = LvalNew;
  LvalHis(iter) = Lval;
  if (valConverge > iterationThres)
    % disp('Finish: Converge on Validation Set.');
    break;
  end
  
  seg = 1;
  startIdx = 1;
  stopIdx  = 10000;
  DBRval = 0;
  while(stopIdx<valLen)
    DBRval   = DBRval + dbr(valLambdaYpre(startIdx:stopIdx)', valY(startIdx:stopIdx)');
    seg      = seg + 1;
    startIdx = 5000*(seg-1)+1;
    stopIdx  = 5000*(seg+1);
  end
  DBRval = DBRval/(seg-1);
  
  if (verbose <= 1)
    plotData(valY(1:10000), valYpre(1:10000), valLambdaYpre(1:10000), LvalHis(1:iter), W)
  end
end

if (verbose <= 2)
  disp([num2str(iter-1),'/',num2str(maxIterations), ...
    ':Train Complete...Lval:',num2str(Lval), 9, '...H:', num2str(H), 9, ...
    '...alpha:',num2str(alpha), 9, '...mu=',num2str(mu)]);
end
W = Whis(bestIter,:);
%% Test data
testLambdaYpre = ANNmodel(testX, W, H, Nx, Nz);
L = logLikelyhood(testY, testLambdaYpre, alpha*norm(W, 1));

seg = 1;
startIdx = 1;
stopIdx  = 10000;
DBR = 0;
while(stopIdx<testLen)
  DBR      = DBR + dbr(testLambdaYpre(startIdx:stopIdx)', testY(startIdx:stopIdx)');
  seg      = seg + 1;
  startIdx = 5000*(seg-1)+1;
  stopIdx  = 5000*(seg+1);
end
DBR = DBR/(seg-1);
if (verbose <= 2)
  disp(['      Test Complete...L:',num2str(L), 9, '...DBR:',num2str(DBR)]);
end
end
