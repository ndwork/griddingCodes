

function testPowerIteration
  clear; close all; rng(1);

  N = 3;
  M = rand(N,N);

  disp(['Norm of M is: ', num2str(norm(M,2))]);
  
  maxIters = 1000;
  tolerance = 0;
  piNorm = powerIteration_dan( M, rand(3,1), maxIters, tolerance );
  
  disp(['Power Iteration norm: ', num2str(piNorm)]);
end


function [nrm,lambda,flag] = powerIteration( M, x, maxIters, tolerance )

  if numel(x) == 0
    sM = size(M);
    x = rand(sM(2),1);
  end

  lambdaPrev = 0;
  lambdaVals = zeros(maxIters,1);
  for iter = 1:maxIters
    MtMx = M'*M*x;
    lambda = norm(MtMx,2);
    if lambda==0, break; end;

    x = MtMx / lambda;

    lambdaVals(iter) = lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      flag = 0;
      break;
    end
  end

  lambdaVals = lambdaVals(1:iter);
  nrm = sqrt( lambda );
end



function [nrm,lambdaVals] = powerIteration_dan(A, x, maxIters, tolerance)
  % Daniel's code

  if numel(x) == 0
    sM = size(M);
    x = rand(sM(2),1);
  end
  x = x/norm(x);

  lambdaPrev = 0;
  lambdaVals = [];

  for iter = 1:maxIters

     Ax = A*x;
     lambda = norm(Ax);
     x = Ax/lambda;

     diff = abs(lambda - lambdaPrev);
     lambdaPrev = lambda;  
     lambdaVals = [lambdaVals,lambda];

  end

  nrm = sqrt(lambda);
end

