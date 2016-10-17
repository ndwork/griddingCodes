function [weights,flag,res] = makePrecompWeights_3D( ...
  traj, N, varargin )
  % [weights,flag,res] = makePrecompWeights_3D( traj, N, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC, 'alg', alg ] )
  %
  % Determine the density pre-compensation weights to be used in gridding
  %
  % Inputs:
  %   traj is a Mx3 element array specifying the k-space trajectory.
  %     The first/second column is the ky/kx/kz location.
  %     The units are normalized to [-0.5,0.5)
  %   N is a 2 element array [Ny Nx Nz] representing the number of grid points
  %
  % Optional Inputs:
  %   alpha - the oversampling factor > 1
  %   W - the window width in pixels
  %   nC - the number of points to sample the convolution kernel
  %   alg - a string specifying the algorithm to use
  %     ls (default) - specifies least squares
  %     fp - specifies fixed point iteration
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
  defaultAlg = 'ls';
  defaultNIter = 15;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
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
    case 'cls'
      % Least squares in frequency domain on grid points
      [weights,flag,res] = makePrecompWeights_3D_CLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );
      
    case 'fdLSDC'
      [weights,flag,res] = makePrecompWeights_3D_fdLSDC( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'fdLSDC_con'
      [weights,flag,res] = makePrecompWeights_3D_fdLSDC_con( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );
      
    case 'fdLSDC_hub'
      [weights,flag,res] = makePrecompWeights_3D_fdLSDC_hub( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'fp'
      % Pipe's method
      [weights,flag,res] = makePrecompWeights_3D_FP( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC, 'nIter', nIter );

    case 'mCLS'
      % CLS masked over region of support in frequency domain
      [weights,flag,res] = makePrecompWeights_3D_mCLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'sdLSDC'
      % Least squares in the space domain.
      [weights,flag,res] = makePrecompWeights_3D_sdLSDC( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'sdLSDC_con'
      % Constrained least squares in the space domain.
      [weights,flag,res] = makePrecompWeights_3D_sdLSDC_con( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'sdLSDC_hub'
      % Constrained least squares in the space domain.
      [weights,flag,res] = makePrecompWeights_3D_sdLSDC_hub( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    otherwise
      error('makePrecompWeights: Algorithm not recognized');
  end

end



function [weights,lsFlag,lsRes] = makePrecompWeights_3D_CLS( ...
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
      out = applyCT_3D( in, traj, cGrid, kCy, kCx, Cy, Cx );
    else
      out = applyC_3D( in, traj, cGrid, kCy, kCx, Cy, Cx );
      out = out(:);

      disp(['makePrecompWeights_3D_CLS working on iteration ', num2str(iteration) ]);
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


function [weights,flag,res] = makePrecompWeights_3D_FP( ...
  traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultNIter = 20;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
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
  Ny=N(1);  Nx=N(2);  Nz=N(3);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );
  Gz = Nz;
  [kCz,Cz,~] = makeKbKernel( Gz, Nz, trueAlpha, W, nC );

  nTraj = size( traj, 1 );
  weights = ones( nTraj, 1 );

  flag = 1;
  for iteration=1:nIter
    if mod( iteration, 5 ) == 0,
      disp(['makePrecompWeights_3D_FP working on iteration ', ...
        num2str(iteration) ]);
    end

    oldWeights = weights;
    denom = applyC_3D( oldWeights, traj, N, kCy, kCx, kCz, Cy, Cx, Cz, traj );
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


function [weights,lsFlag,lsRes] = makePrecompWeights_3D_mCLS( ...
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

  cGrid = 2*N;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, cGrid );
      out = applyCT_3D( in, traj, cGrid, kCy, kCx, Cy, Cx );
    else
      out = applyC_3D( in, traj, cGrid, kCy, kCx, Cy, Cx );
      out = out(:);

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_3D_mCLS working on iteration ', ...
          num2str(iteration) ]);
      end
      
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b=ones(cGrid);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N, 'imgTitle', 'mCLS');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end



function [weights,lsFlag,lsRes] = makePrecompWeights_3D_fdLSDC( ...
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
      out = applyCT_3D( in, traj, 0, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_3D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_3D_fdLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b = ones(size(traj,1),1);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), ...
    tolerance, maxIter );

  scale = showPSF( weights, traj, N, 'imgTitle', 'fdLSDC');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end


function [weights,flag,residual] = makePrecompWeights_3D_fdLSDC_con( ...
  traj, N, varargin )
  % Least squares analogy of Pipe's algorithm

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
      out = applyCT_3D( in, traj, 0, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_3D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_3D_fdLSDC_con working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  b = ones(size(traj,1),1);
  nTraj = size(traj,1);
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
  opts.maxIts = 1000;
  opts.printEvery = 1;
  opts.tol = 1d-4;
  x0 = ones( nTraj, 1 );
  weights = tfocs( smooth_quad, { @linop_4_tfocs, -b }, proj_Rplus, x0, opts );
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




function [weights,flag,residual] = makePrecompWeights_3D_fdLSDC_hub( ...
  traj, N, varargin )
  % Least squares analogy of Pipe's algorithm

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
      out = applyCT_3D( in, traj, 0, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_3D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_3D_fdLSDC_hub working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  b = ones(size(traj,1),1);
  nTraj = size(traj,1);
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
  opts.maxIts = 1000;
  opts.printEvery = 1;
  opts.tol = 1d-4;
  x0 = ones( nTraj, 1 );
  weights = tfocs( smooth_huber(0.02), { @linop_4_tfocs, -b }, [], x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end
  if nargout > 2
    psf = applyA( weights );
    residual = norm( psf(:) - b(:), 2 ) / norm(b(:),2);
  end
  
tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  scale = showPSF( weights, traj, N, 'imgTitle', 'fdLSDC_hub');
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N );
end



function [weights,lsFlag,lsRes] = makePrecompWeights_3D_sdLSDC( ...
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

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  cgrid = 2*N;
  radialImg = makeRadialImg( cgrid );
  b=zeros(cgrid);  b(1,1)=1;  b=fftshift(b);
  mask = double( radialImg <= min(N) );
  normWeights = ones(cgrid);
  normWeights(b>0) = sqrt( prod(cgrid) );
  b = normWeights .* b;

  iteration = 0;
  function out = applyWA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      nIn = numel(in);
      in = in(1:nIn/2) + 1i*(in(nIn/2+1:end));
      in = reshape( in, cgrid );
      masked = mask .* normWeights .* in;
      out = iGrid_3D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      iGridTed = iGridT_3D( in, traj, cgrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      masked = normWeights .* mask .* iGridTed;
      out = [ real(masked(:)); imag(masked(:)); ];

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0
        disp(['makePrecompWeights_3D_sdLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  tolerance = 1d-6;
  maxIter = 1000;
  Wb = [ real(b(:)); imag(b(:)); ];
  [weights,lsFlag,lsRes] = lsqr( @applyWA, Wb(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N, mask, 'imgTitle', 'sdLSDC' );
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end


function [weights,flag,residual] = makePrecompWeights_3D_sdLSDC_con( ...
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

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  cgrid = 2*N;
  radialImg = makeRadialImg( cgrid );
  b=zeros(cgrid);  b(1,1)=1;  b=fftshift(b);
  mask = double( radialImg <= min(N) );
  normWeights = ones(cgrid);
  normWeights(b>0) = sqrt( prod(cgrid) );
  b = normWeights .* b;

  iteration = 0;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      nIn = numel(in);
      in = in(1:nIn/2) + 1i*(in(nIn/2+1:end));
      in = reshape( in, cgrid );
      masked = mask .* normWeights .* in;
      out = iGrid_3D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      iGridTed = iGridT_3D( in, traj, cgrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      masked = normWeights .* mask .* iGridTed;
      out = [ real(masked(:)); imag(masked(:)) ];

      iteration = iteration + 1;
      if mod( iteration, 5 )==0
        disp(['makePrecompWeights_3D_con_sdLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  nTraj = size(traj,1);
  Wb = [ real(b(:)); imag(b(:)); ];
  function y = linop_4_tfocs( in, mode )
    switch mode,
      case 0
        y = [numel(Wb), nTraj];
      case 1
        y = applyA( in );
      case 2
        y = applyA( in, 'transp' );
    end
  end
  %varargout = linop_test( @linop_4_tfocs );

  opts.alg = 'N83';
  opts = tfocs;
  opts.maxIts = 1000;
  opts.printEvery = 1;
  opts.tol = 1d-4;
  x0 = ones( nTraj, 1 );
  %[ x, out ] = tfocs( smoothF, affineF, nonsmoothF, x0, opts );
  %weights = tfocs_N83( smooth_quad, { @linop_4_tfocs, -Wb(:) }, proj_Rplus, x0, opts );
  %weights = tfocs( smooth_quad, { @linop_4_tfocs, -Wb(:) }, proj_Rplus, x0, opts );
  weights = tfocs( smooth_huber(0.02), { @linop_4_tfocs, -Wb(:) }, proj_Rplus, x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end
  if nargout > 2
    Wpsf = applyA( weights );
    residual = norm( Wpsf(:) - Wb(:), 2 ) / norm(Wb(:),2);
  end

  scale = showPSF( weights, traj, N, mask );
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end



function [weights,flag,residual] = makePrecompWeights_3D_sdLSDC_hub( ...
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

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  cgrid = 2*N;
  radialImg = makeRadialImg( cgrid );
  b=zeros(cgrid);  b(1,1)=1;  b=fftshift(b);
  mask = double( radialImg <= min(N) );
  normWeights = ones(cgrid);
  normWeights(b>0) = sqrt( prod(cgrid) );
  b = normWeights .* b;

  iteration = 0;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      nIn = numel(in);
      in = in(1:nIn/2) + 1i*(in(nIn/2+1:end));
      in = reshape( in, cgrid );
      masked = mask .* normWeights .* in;
      out = iGrid_3D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      iGridTed = iGridT_3D( in, traj, cgrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      masked = normWeights .* mask .* iGridTed;
      out = [ real(masked(:)); imag(masked(:)) ];

      iteration = iteration + 1;
      if mod( iteration, 5 )==0
        disp(['makePrecompWeights_3D_con_sdLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  nTraj = size(traj,1);
  Wb = [ real(b(:)); imag(b(:)); ];
  function y = linop_4_tfocs( in, mode )
    switch mode,
      case 0
        y = [numel(Wb), nTraj];
      case 1
        y = applyA( in );
      case 2
        y = applyA( in, 'transp' );
    end
  end
  %varargout = linop_test( @linop_4_tfocs );

  opts.alg = 'N83';
  opts = tfocs;
  opts.maxIts = 1000;
  opts.printEvery = 1;
  opts.tol = 1d-4;
  x0 = ones( nTraj, 1 );
  %[ x, out ] = tfocs( smoothF, affineF, nonsmoothF, x0, opts );
  %weights = tfocs_N83( smooth_quad, { @linop_4_tfocs, -Wb(:) }, proj_Rplus, x0, opts );
  %weights = tfocs( smooth_quad, { @linop_4_tfocs, -Wb(:) }, proj_Rplus, x0, opts );
  weights = tfocs( smooth_huber(0.02), { @linop_4_tfocs, Wb(:) }, [], x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end
  if nargout > 2
    Wpsf = applyA( weights );
    residual = norm( Wpsf(:) - Wb(:), 2 ) / norm(Wb(:),2);
  end

  scale = showPSF( weights, traj, N, mask );
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end


