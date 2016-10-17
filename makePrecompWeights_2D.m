
function [weights,flag,res] = makePrecompWeights_2D( ...
  traj, N, varargin )
  % [weights,flag,res] = makePrecompWeights_2D( traj, N, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC, 'alg', alg ] )
  %
  % Determine the density pre-compensation weights to be used in gridding
  %
  % Inputs:
  %   traj is a Mx2 element array specifying the k-space trajectory.
  %     The first/second column is the ky/kx location.
  %     The units are normalized to [-0.5,0.5)
  %   N is a 2 element array [Ny Nx] representing the number of grid points
  %
  % Optional Inputs:
  %   alpha - the oversampling factor > 1
  %   W - the window width in pixels
  %   nC - the number of points to sample the convolution kernel
  %   alg - a string specifying the algorithm to use
  %     rLSDC (default) - robust least squares
  %     LSDC - specifies least squares
  %     FP - specifies fixed point iteration
  %     CLS - minimizes in the frequency domain
  %   nIter - specifies the number of iterations of fp method
  %
  % Outputs:
  %   weights - 1D array with density compensation weights
  %
  % Optional Outputs:
  %   flag - flag describing results of optimization (see lsqr
  %     documentation)
  %   res - residual
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultAlg = 'rLSDC';
  defaultNIter = 15;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.addParameter( 'alg', defaultAlg );
  p.addParameter( 'nIter', defaultNIter, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  alg = p.Results.alg;
  nIter = p.Results.nIter;

  flag = 0;
  res = 0;
  switch alg
    case 'CLS'
      % Least squares in frequency domain on grid points
      [weights,flag,res] = makePrecompWeights_2D_CLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'FP'
      % Pipe's method
      [weights,flag,res] = makePrecompWeights_2D_FP( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC, 'nIter', nIter );
      
    case 'LSDC'
      % Optimization analogy of Pipe's Method
      [weights,flag,res] = makePrecompWeights_2D_LSDC( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'rLSDC'
      % LSDC with non-negativity constraint
      [weights,flag,res] = makePrecompWeights_2D_rLSDC( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    otherwise
      error('makePrecompWeights: Algorithm not recognized');
  end

end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_CLS( ...
  traj, N, varargin )

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  iteration = 0;

  % Make the Kaiser Bessel convolution kernel
  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  %cGrid = 2*N;
  cGrid = nGrid;

  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, cGrid );
      out = applyCT_2D( in, traj, cGrid, kCy, kCx, Cy, Cx );
    else
      out = applyC_2D( in, traj, cGrid, kCy, kCx, Cy, Cx );
      out = out(:);

      disp(['makePrecompWeights_2D_CLS working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b = ones(cGrid);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), ...
    tolerance, maxIter );

  radialImg = makeRadialImg( 2*N );
  mask = double( radialImg <= min(N) );
  scale = showPSF( weights, traj, N, mask, 'imgTitle', 'cls');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end


function [weights,flag,res] = makePrecompWeights_2D_FP( ...
  traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultNIter = 20;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.addParameter( 'nIter', defaultNIter, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  nIter = p.Results.nIter;

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  % Make the Kaiser Bessel convolution kernel
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  nTraj = size( traj, 1 );
  weights = ones( nTraj, 1 );

  flag = 1;
  for iteration=1:nIter
    if mod( iteration, 5 ) == 0,
      disp(['makePrecompWeights_2D_FP working on iteration ', ...
        num2str(iteration) ]);
    end

    oldWeights = weights;
    denom = applyC_2D( oldWeights, traj, N, kCy, kCx, Cy, Cx, traj );
    weights = oldWeights ./ denom;
  end

  radialImg = makeRadialImg( 2*N );
  mask = double( radialImg <= min(N) );

  scale = showPSF( weights, traj, N, mask, 'imgTitle', 'fp' );
  weights = scale * weights;

  if nargout > 1
    flag = 0;
  end
  if nargout > 2
    res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );
  end
end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_LSDC( ...
  traj, N, varargin )
  % Optimization analogy of Pipe's algorithm

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  iteration = 0;

  % Make the Kaiser Bessel convolution kernel
  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  cGrid = 2*N;
  
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      out = applyCT_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_2D_LSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b = ones(size(traj,1),1);
  tolerance = 1d-5;
  maxIter = 5000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), ...
    tolerance, maxIter );

  scale = showPSF( weights, traj, N, 'imgTitle', 'fdLSDC');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end


function [weights,flag,residual] = makePrecompWeights_2D_rLSDC( ...
  traj, N, varargin )
  % Method of LSDC with a non-negativity constraint

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  iteration = 0;

  % Make the Kaiser Bessel convolution kernel
  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      out = applyCT_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_2D_rLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  nTraj = size(traj,1);
  b = ones(nTraj,1);
  function y = linop_4_tfocs( in, mode )
    switch mode,
      case 0
        y = [numel(b), nTraj];
      case 1
        y = applyA( in );
      case 2
        y = applyA( in, 'transp' );
    end
  end
  %varargout = linop_test( @linop_4_tfocs );

  opts.alg = 'N83';
  opts = tfocs;
  opts.maxIts = 10000;
  opts.printEvery = 1;
  opts.tol = 1d-5;
  x0 = ones( nTraj, 1 );
  [weights,tfocsOptOut] = tfocs( smooth_quad, { @linop_4_tfocs, -b }, ...
    proj_Rplus, x0, opts );
  %weights = tfocs( smooth_huber(0.02), { @linop_4_tfocs, -b }, [], x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end
  if nargout > 2
    psf = applyA( weights );
    residual = norm( psf(:) - b(:), 2 ) / norm(b(:),2);
  end
  
tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  scale = showPSF( weights, traj, N, 'imgTitle', 'fdLSDC_con');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end

