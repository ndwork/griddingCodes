
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
    case 'cls'
      [weights,flag,res] = makePrecompWeights_2D_CLS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'fp'
      [weights,flag,res] = makePrecompWeights_2D_FP( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'ls'
      [weights,flag,res] = makePrecompWeights_2D_LS( ...
        traj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'rls'
      % Robust least squares
      [weights,flag,res] = makePrecompWeights_2D_RLS( ...
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
  kws = [ kwy kwx ];

  %nGrid = 2*N;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, nGrid );
      out = applyCT_2D( in, traj, nGrid, kws, kCy, kCx, Cy, Cx );
    else
      out = applyC_2D( in, traj, nGrid, kws, kCy, kCx, Cy, Cx );
      out = out(:);

      disp(['makePrecompWeights_2D working on iteration ', num2str(iteration) ]);
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

psf = iGridT_2D( weights, traj, 2*N, 'alpha', 2.0, 'W', W, 'nC', nC );
radialImg = makeRadialImg( size(psf) );
mask = radialImg <= min(N);
psf = mask .* psf;
psf = cropData( psf, 2*N );
psf = psf ./ max(psf(:));
b=zeros(size(psf)); b(1,1)=1; b=fftshift(b);
mse = sum( abs( mask(:).*psf(:) - mask(:).*b(:) ).^2 ) ./ sum(mask(:));
disp(['MSE: ', num2str(mse)]);
tmp=abs(psf);  tmp=tmp./max(tmp(:)); tmp=20*log10(tmp);
tmp(~isfinite(tmp))=0; imshow( imresize(tmp,3,'nearest'), [-100 0] );
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
  [kCy,Cy,~,kwy] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~,kwx] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );
  kws = [ kwy kwx ];

  nTraj = size( traj, 1 );
  weights = ones( nTraj, 1 );
  maxIteration = 20;
  tolerance = 1d-5;
  flag = 1;
  for iteration=1:maxIteration
    disp(['makePrecompWeights_2D working on iteration ', num2str(iteration) ]);

    oldWeights = weights;
    denom = applyC_2D( oldWeights, traj, N, kws, kCy, kCx, Cy, Cx, traj );
    weights = oldWeights ./ denom;

    res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );
    if res < tolerance
      flag = 0;   % the fixed point iteration has converged
      break;
    end
  end

psf = iGridT_2D( weights, traj, 2*N, 'alpha', 2.0, 'W', W, 'nC', nC );
radialImg = makeRadialImg( size(psf) );
mask = radialImg <= min(N);
psf = mask .* psf;
psf = cropData( psf, 2*N );
psf = psf ./ max(psf(:));
b=zeros(size(psf)); b(1,1)=1; b=fftshift(b);
mse = sum( abs( mask(:).*psf(:) - mask(:).*b(:) ).^2 ) ./ sum(mask(:));
disp(['MSE: ', num2str(mse)]);
tmp=abs(psf);  tmp=tmp./max(tmp(:)); tmp=20*log10(tmp);
tmp(~isfinite(tmp))=0; imshow( imresize(tmp,3,'nearest'), [-100 0] );

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

  nGrid = 2.5 * N;
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

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes] = lsqr( @applyA, b(:), tolerance, maxIter );

psf = iGridT_2D( weights, traj, nGrid, 'alpha', 2.0, 'W', W, 'nC', nC );
psf = cropData( psf, 2*N );
mask = cropData( mask, 2*N );
psf = mask .* psf;
psf = psf ./ max(psf(:));
b=zeros(size(psf)); b(1,1)=1; b=fftshift(b);
mse = sum( abs( mask(:).*psf(:) - mask(:).*b(:) ).^2 ) ./ sum(mask(:));
disp(['MSE: ', num2str(mse)]);
tmp=abs(psf);  tmp=tmp./max(tmp(:)); tmp=20*log10(tmp);
tmp(~isfinite(tmp))=0; imshow( imresize(tmp,3,'nearest'), [-100 0] );
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


