
function weights = makePrecompWeights_1D( traj, N, varargin )
  % weights = makePrecompWeights_1D( traj, N, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Determine the density pre-compensation weights to be used in gridding
  %
  % Inputs:
  %   traj is a M element array specifying the k-space trajectory.
  %     The units are normalized to [-0.5,0.5).
  %   N is the number of grid points
  %
  % Optional Inputs:
  %   alpha is the oversampling factor > 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Written by Nicholas Dwork - Copyright 2016

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
  trueAlpha = nGrid / N;

  % Make the Kaiser Bessel convolution kernel
  G = nGrid;
  [kC,C,c1D,kw] = makeKbKernel( G, nGrid, trueAlpha, W, nC );

  iteration = 0;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      iteration = iteration + 1;
      if mod( iteration, 10 ) == 0, disp(['lsqr iteration: ', num2str(iteration)]); end;
      %out = applyCT_1D( in, traj, nGrid, kw, kC, C );
      out = 1/nGrid * iGrid_1D( in, traj, 'alpha', alpha, 'W', W, 'nC', nC );
    else
      %out = applyC_1D( in, traj, nGrid, kw, kC, C );
      out = 1/nGrid * iGridT_1D( in, traj, nGrid, 'alpha', alpha, 'W', W, 'nC', nC );
    end
  end

  %ks = size2fftCoordinates( nGrid );
  %b = applyC_1D( ones(nGrid,1), ks, nGrid, kw, kC, C );
  %b = zeros(nGrid,1);  b(1) = 1;  b=fftshift(b);
  %b = 1/trueAlpha * ones(nGrid,1);
  b = zeros(nGrid,1);  b(1) = 1/trueAlpha;  b=fftshift(b);
  weights = lsqr( @applyA, b, 1d-5, 1000 );
end

