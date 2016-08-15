function [weights,flag,res] = makePrecompWeights_2D_v2( ...
  traj, N, varargin )
  % [weights,flag,res] = makePrecompWeights_2D( traj, N, ...
  %   [ ros, 'alpha', alpha, 'W', W, 'nC', nC, 'alg', alg ] )
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
  %   ros - the Region of Support
  %   alpha - the oversampling factor > 1
  %   W - the window width in pixels
  %   nC - the number of points to sample the convolution kernel
  %   alg - a string specifying the algorithm to use
  %     cls - minimizes in the frequency domain
  %     ls (default) - specifies least squares
  %     rls - robust least squares
  %     fp - specifies fixed point iteration
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

  defaultROS = [];
  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultAlg = 'ls';
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addOptional( 'ros', defaultROS );
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.addParamValue( 'alg', defaultAlg );
  p.parse( varargin{:} );
  ros = p.Results.ros;
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  alg = p.Results.alg;

  flag = 0;
  res = 0;
  switch alg
    case 'cls'
      [weights,flag,res] = makePrecompWeights_2D_CLS( ...
        traj, N, ros, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'fp'
      [weights,flag,res] = makePrecompWeights_2D_FP( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'mCLS'
      % CLS masked over region of support in frequency domain
      [weights,flag,res] = makePrecompWeights_2D_mCLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'regCLS'
      [weights,flag,res] = makePrecompWeights_2D_regCLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'regWLS'
      % Regularized least squares determination of the weights in the space domain.
      [weights,flag,res] = makePrecompWeights_2D_regWLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'wls'
      % Least squares determination of the weights in the space domain.
      [weights,flag,res] = makePrecompWeights_2D_WLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );
      
    otherwise
      error('makePrecompWeights: Algorithm not recognized');
  end

end



function [weights,lsFlag,lsRes] = makePrecompWeights_2D_CLS( ...
  traj, N, varargin )

  defaultROS = [];
  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addOptional( 'ros', defaultROS );
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  ros = p.Results.ros;
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
  
  if isempty( ros )
    ros = ones( cGrid );
  end
  
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, cGrid );
      in = in .* ros;
      out = applyCT_2D( in, traj, cGrid, kCy, kCx, Cy, Cx );
    else
      out = applyC_2D( in, traj, cGrid, kCy, kCx, Cy, Cx );
      out = out .* ros;
      out = out(:);

      disp(['makePrecompWeights_2D_CLS working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b = ones(cGrid) .* ros;
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), ...
    tolerance, maxIter );

  scale = showPSF( weights, traj, N, 'imgTitle', 'cls');
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
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

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
  maxIteration = 20;


  tolerance = 1d-5;
  flag = 1;
  for iteration=1:maxIteration
    disp(['makePrecompWeights_2D_FP working on iteration ', num2str(iteration) ]);

    oldWeights = weights;
    denom = applyC_2D( oldWeights, traj, N, kCy, kCx, Cy, Cx, traj );
    weights = oldWeights ./ denom;

    res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );
    if res < tolerance
      flag = 0;   % the fixed point iteration has converged
      break;
    end
  end

  scale = showPSF( weights, traj, N, 'imgTitle', 'fp' );
  weights = scale * weights;
end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_mCLS( ...
  traj, N, varargin )

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
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

  b=ones(cGrid);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N, 'imgTitle', 'cls');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end



function [weights,lsFlag,lsRes] = makePrecompWeights_2D_regCLS( ...
  traj, N, varargin )

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  sqrtGamma = sqrt(1);
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
      in1 = in(1:prod(cGrid));
      tmp1 = applyCT_2D( in1, traj, cGrid, kCy, kCx, Cy, Cx );

      in2 = in(prod(cGrid)+1:end);
      tmp2 = sqrtGamma * in2;
      
      out = tmp1 + tmp2;
    else
      out1 = applyC_2D( in, traj, cGrid, kCy, kCx, Cy, Cx );
      out2 = sqrtGamma * in;
      out = [ out1(:); out2(:); ];

      disp(['makePrecompWeights_2D_CLS working on iteration ', ...
        num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp1=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  nTraj = size( traj, 1 );
  b1=ones(cGrid);
  b2 = zeros( nTraj, 1 );
  b = [ b1(:); b2(:); ];
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N , 'imgTitle', 'regCLS');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_regWLS( ...
  traj, N, varargin )

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  cgrid = 2*N;
  sqrtGamma = sqrt(1);
  iteration = 0;
  function out = applyWA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in1 = in(1:prod(cgrid));
      in1 = reshape( in1, cgrid );
      masked = mask .* normWeights1 .* in1;
      tmp1 = iGrid_2D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      tmp1 = real(tmp1);

      in2 = in(prod(cgrid)+1:end);
      tmp2 = sqrtGamma * in2;

      out = tmp1 + tmp2;
    else
      out1 = iGridT_2D( in, traj, cgrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out1 = normWeights1 .* mask .* out1;
      
      out2 = sqrtGamma * in;

      out = [ out1(:); out2(:); ];

      disp(['makePrecompWeights_2D_regWLS working on iteration ', ...
        num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  radialImg = makeRadialImg( cgrid );
  b1=zeros(cgrid);  b1(1,1)=1;  b1=fftshift(b1);
  mask = double( radialImg <= min(N) );
  normWeights1 = ones(cgrid);
  normWeights1(b1>0) = sqrt( prod(cgrid) );
  b1 = normWeights1 .* b1;
  nTraj = size( traj, 1 );
  b2 = zeros( nTraj, 1 );
  b = [ b1(:); b2(:); ];
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes] = lsqr( @applyWA, b(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N, mask );
oldWeights = weights;
  weights = scale * weights;
scale = showPSF( weights, traj, N, mask );

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end



function [weights,lsFlag,lsRes] = makePrecompWeights_2D_WLS( ...
  traj, N, varargin )

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  cgrid = 2*N;
  radialImg = makeRadialImg( cgrid );
  b=zeros(cgrid);  b(1,1)=1;  b=fftshift(b);
  mask = double( radialImg <= min(N) );
  mask(b>0) = sqrt( prod(cgrid) );
  b(b>0) = sqrt( prod(cgrid) );

  iteration = 0;
  function out = applyWA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, cgrid );
      masked = mask .* in;
      out = iGrid_2D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      out = iGridT_2D( in, traj, cgrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = mask .* out;
      out = out(:);

      disp(['makePrecompWeights_2D_WLS working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes] = lsqr( @applyWA, b(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N, mask );
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end



