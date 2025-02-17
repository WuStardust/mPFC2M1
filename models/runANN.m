function [W,L,DBR,Lval,LtrainHis,LvalHis,DBRtrainHis,DBRvalHis,Whis,muHis,HessianDetHis] = runANN(H, Nz, xi1, xi2, mu, threshold, ...
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
DBRvalHis = zeros(1, maxIterations);
DBRtrainHis = zeros(1, maxIterations);
HessianDetHis = zeros(1, maxIterations);
muHis = zeros(1, maxIterations);
%% train ANN
% forward
[trainLambdaYpre, ~, trainLambdaZpre] = ANNmodel(trainX, W, Nx, Nz);
Ltrain = logLikelyhood(trainY, trainLambdaYpre, alpha*norm(W, 1));
LtrainHis(iter) = Ltrain;
valLambdaYpre= ANNmodel(valX, W, Nx, Nz);
Lval = logLikelyhood(valY, valLambdaYpre, alpha*norm(W, 1));
LvalHis(iter) = Lval;
LvalBest = Lval;
bestIter = iter;

DBRtrain = adbr(trainLambdaYpre, trainY, trainLen);
DBRval = adbr(valLambdaYpre, valY, valLen);

while (iter<maxIterations)
  DBRtrainHis(iter) = DBRtrain;
  DBRvalHis(iter) = DBRval;
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
  HessianDetHis(iter) = det(Hessian);
  if (rcond(Hessian) < 1e-15 || sum(sum(isnan(Hessian))) == 1)
    disp('Error: BAD He!');
    break;
  end
  
  while(1)
    % lm update
    Wnew = W - Grad/(Hessian-mu*eye(size(Hessian)));

    % forward
    [trainLambdaYpre, ~, trainLambdaZpre] = ANNmodel(trainX, Wnew, Nx, Nz);
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
  muHis(iter) = mu;
  
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
  [valLambdaYpre, valYpre] = ANNmodel(valX, W, Nx, Nz);
  LvalNew = logLikelyhood(valY, valLambdaYpre, alpha*norm(W, 1));
  if (LvalBest <= LvalNew)
    LvalBest = LvalNew;
    bestIter = iter;
  end
  LvalBest = max(LvalBest, LvalNew);
  % L on validation set change too little, or drop too much
  if (abs(LvalNew-Lval)<threshold || LvalNew-LvalBest<-500)
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
  
  DBRtrain = adbr(trainLambdaYpre, trainY, trainLen);
  DBRval = adbr(valLambdaYpre, valY, valLen);

  if (verbose <= 1)
    plotData(valY, valYpre, valLambdaYpre, LvalHis(1:iter), W)
  end
end

if (verbose <= 2)
  disp([num2str(iter-1),'/',num2str(maxIterations), ...
    ':Train Complete...Lval:',num2str(Lval), 9, '...H:', num2str(H), 9, ...
    '...Nz:', num2str(Nz), 9, '...alpha:',num2str(alpha), 9, '...mu=',num2str(mu)]);
end
W = Whis(bestIter,:);
LtrainHis(LtrainHis==0) = nan;
LvalHis(LvalHis==0) = nan;
Whis = Whis(1:iter,:);
%% Test data
testLambdaYpre = ANNmodel(testX, W, Nx, Nz);
L = logLikelyhood(testY, testLambdaYpre, alpha*norm(W, 1));

DBR = adbr(testLambdaYpre, testY, testLen);

if (verbose <= 3)
  disp(['Test Complete...L:',num2str(L), 9, '...H:', num2str(H), 9, ...
    '...Nz:', num2str(Nz), 9, '...xi1:', num2str(xi1), 9, '...xi2:', num2str(xi2), 9, ...
    '...DBR:',num2str(DBR)]);
end
end
