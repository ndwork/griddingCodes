
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

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultAlg = 'ls';
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.addParamValue( 'alg', defaultAlg );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  alg = p.Results.alg;

  flag = 0;
  res = 0;
  switch alg
    case '2cls'
      [weights,flag,res] = makePrecompWeights_2D_2CLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'cls'
      [weights,flag,res] = makePrecompWeights_2D_CLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'fp'
      [weights,flag,res] = makePrecompWeights_2D_FP( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'ls'
      [weights,flag,res] = makePrecompWeights_2D_WLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'rls_cp'
      % Robust least squares
      [weights,flag,res] = makePrecompWeights_2D_RLSCP( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'rls'
      % Robust least squares with FISTA
      [weights,flag,res] = makePrecompWeights_2D_WRLSF( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );
      
    case 'vor'
      weights = voronoidens(traj);

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
  [kCy,Cy,~,kwy] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~,kwx] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  %nGrid = 2*N;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, nGrid );
      out = applyCT_2D( in, traj, nGrid, kCy, kCx, Cy, Cx );
    else
      out = applyC_2D( in, traj, nGrid, kCy, kCx, Cy, Cx );
      out = out(:);

      disp(['makePrecompWeights_2D_CLS working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b=ones(nGrid);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes] = lsqr( @applyA, b(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N );
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

  scale = showPSF( weights, traj, N );
  weights = scale * weights;
end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_LS( ...
  traj, N, varargin )

%   defaultAlpha = 1.5;
%   defaultW = 8;
%   defaultNc = 500;
%   checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
%   p = inputParser;
%   p.addParamValue( 'alpha', defaultAlpha, checknum );
%   p.addParamValue( 'W', defaultW, checknum );
%   p.addParamValue( 'nC', defaultNc, checknum );
%   p.parse( varargin{:} );
%   alpha = p.Results.alpha;
%   W = p.Results.W;
%   nC = p.Results.nC;


  iteration = 0;

  nGrid = ceil( 2.5 * N );
  trueAlpha = max( nGrid ./ N );
  W = 8;
  nC = 500;

  radialImg = makeRadialImg( nGrid );
  mask = radialImg <= min(N);

  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, nGrid );
      masked = mask .* in;
      out = iGrid_2D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      out = iGridT_2D( in, traj, nGrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = mask .* out;
      out = out(:);

      disp(['makePrecompWeights_2D_LS working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes] = lsqr( @applyA, b(:), tolerance, maxIter );

  scale = showPSF( weights, traj, N, mask );
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end



function [weights,lsFlag,lsRes] = makePrecompWeights_2D_WLS( ...
  traj, N, varargin )

  iteration = 0;

  nGrid = ceil( 2.5 * N );
  trueAlpha = max( nGrid ./ N );
  W = 8;
  nC = 500;

  radialImg = makeRadialImg( nGrid );
  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b);
  mask = double( radialImg <= min(N) );
  mask(b>0) = sqrt( 2*prod(nGrid) );

  function out = applyWA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, nGrid );
      masked = mask .* in;
      out = iGrid_2D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      out = iGridT_2D( in, traj, nGrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = mask .* out;
      out = out(:);

      disp(['makePrecompWeights_2D_LS working on iteration ', num2str(iteration) ]);
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


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_RLSCP( ...
  traj, N, varargin )

  % Uses Chambolle-Pock to minimize Robust-Least Squares optimization
  % problem for determining density compensation weights

  error('UNTESTED - DOESN''T WORK');
  
  maxIterCP = 1000;
  iteration = 0;

  nGrid = ceil( 2.5 * N );
  trueAlpha = max( nGrid ./ N );
  W = 8;
  nC = 500;

  radialImg = makeRadialImg( nGrid );
  mask = radialImg <= min(N);

  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, nGrid );
      masked = mask .* in;
      out = iGrid_2D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      out = iGridT_2D( in, traj, nGrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = mask .* out;
      out = out(:);

      disp(['makePrecompWeights_2D working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b);
  bNonzeroIndx = find( b ~= 0 );

  x0 = rand( size(mask) );
  
  %[nrmA,piFlag] = powerIteration( @applyA, x0 );
nrmA = 0; piFlag=0; load 'nrmA_64x64.mat';
  nrmASq = nrmA*nrmA;
  tau = 1 / nrmA * 0.98;
  sigma = 1 / nrmA * 0.98;

  weights = ones( size(traj,1), 1 );
  wBar = zeros( size(traj,1), 1 );
  y = applyA( weights );
  alpha = 1;
  residuals = zeros( maxIterCP, 1 );  figure;
  for cpIndx=1:maxIterCP
    % Update y
    tmp = y + sigma * applyA( wBar );
    y = tmp / (sigma + 1);
    y(bNonzeroIndx) = y(bNonzeroIndx) - sigma/(sigma+1);

    % Update w
    lastW = weights;
    ATy = applyA( y, 'transp' );
    weights = max( weights - tau*ATy, 0 );

    % Update wBar
    wBar = weights + alpha * ( weights - lastW );

    % To see convergence, calculate residual
    Aw = applyA( weights );
    residuals(cpIndx) = norm( Aw(:) - b(:), 2 );
    if cpIndx>1
      plot( 1:cpIndx, residuals(1:cpIndx), 'LineWidth', 2 );
      drawnow;
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  scale = showPSF( weights, traj, N, mask );
  weights = scale * weights;
end


function [weights,flag,residual] = makePrecompWeights_2D_WRLSF( ...
  traj, N, varargin )
  % Uses FISTA to minimize Robust-Least Squares optimization
  % problem for determining density compensation weights

  nGrid = ceil( 2.5 * N );
  trueAlpha = max( nGrid ./ N );
  W = 8;
  nC = 500;

  radialImg = makeRadialImg( nGrid );
  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b);
  mask = double( radialImg <= min(N) );
  mask(b>0) = sqrt( 2*prod(nGrid) );

  function out = applyWA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      nIn = numel(in);
      in = in(1:nIn/2) + 1i*(in(nIn/2+1:end));
      in = reshape( in, nGrid );
      masked = mask .* in;
      out = iGrid_2D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      iGridTed = iGridT_2D( in, traj, nGrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      masked = mask .* iGridTed;
      out = [ real(masked(:)); imag(masked(:)) ];
    end
  end

  nTraj = size(traj,1);
  Wb = [ real(b(:)); imag(b(:)); ];

  function y = linop_4_tfocs( in, mode )
    switch mode,
      case 0
        y = [numel(Wb), nTraj];
      case 1
        y = applyWA( in );
      case 2
        y = applyWA( in, 'transp' );
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
  weights = tfocs( smooth_quad, { @linop_4_tfocs, -Wb(:) }, proj_Rplus, x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end

  if nargout > 2
    Wpsf = applyWA( weights );
    residual = norm( Wpsf(:) - Wb(:), 2 ) / norm(Wb(:),2);
  end

  scale = showPSF( weights, traj, N, mask );
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end



function [weights,flag,residual] = makePrecompWeights_2D_RLSF( ...
  traj, N, varargin )

  % Uses FISTA to minimize Robust-Least Squares optimization
  % problem for determining density compensation weights

  nGrid = ceil( 2.5 * N );
  trueAlpha = max( nGrid ./ N );
  W = 8;
  nC = 500;

  radialImg = makeRadialImg( nGrid );
  mask = radialImg <= min(N);

  iteration = 0;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      nIn = numel(in);
      in = in(1:nIn/2) + 1i*(in(nIn/2+1:end));
      in = reshape( in, nGrid );
      masked = mask .* in;
      out = iGrid_2D( masked, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      iGridTed = iGridT_2D( in, traj, nGrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      masked = mask .* iGridTed;
      out = [ real(masked(:)); imag(masked(:)) ];

      disp(['makePrecompWeights_2D_RLS working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  nTraj = size(traj,1);
  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b); b=[real(b(:)); imag(b(:));];

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
  %[ x, out ] = tfocs( smoothF, affineF, nonsmoothF, x0, opts );
  x0 = ones( nTraj, 1 );
  %weights = tfocs_N83( smooth_quad, { @linop_4_tfocs, -b(:) }, proj_Rplus, x0, opts );
  weights = tfocs( smooth_quad, { @linop_4_tfocs, -b(:) }, proj_Rplus, x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end

  if nargout > 2
    psf = applyA( weights );
    residual = norm( psf(:) - b(:), 2 ) / norm(b(:),2);
  end

  scale = showPSF( weights, traj, N, mask );
  weights = scale * weights;

  %noise = randn(size(weights))*0.1;
  %noisyWeights = weights .* ( 1 + noise );
  %showPSF( noisyWeights, traj, N, mask );
end



function areas = voronoidens(kTraj)
  %
  % function areas = voronoidens(kTraj);
  %
  % output: areas of cells for each point 
  %           (if point doesn't have neighbors the area is NaN)
  %
  % Written by John Pauly

  [V,C] = voronoin(kTraj);  % returns vertices and cells of voronoi diagram
  nK = size( kTraj, 1 );
  areas = zeros( nK, 1 );
  for j = 1:nK
    x = V(C{j},1); y = V(C{j},2); lxy = length(x);
    A = abs( sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:)) ) );
    areas(j) = A;
  end

  areas( ~isfinite(areas) ) = 0;

  %ky=kTraj(:,1);  kx=kTraj(:,2);
  %[vx, vy] = voronoi(kx,ky);
  %figure; plot(kx,ky,'r.',vx,vy,'b-'); axis equal
end


