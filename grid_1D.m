
function out = grid_1D( F, traj, N, varargin )
  % F = grid_1D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI encoding with Inverse Gridding
  %
  % Inputs
  %   F is a 1D array representing the Fourier values
  %   traj is a M element array specifying the k-space trajectory.
  %     The units are normalized to [-0.5,0.5).
  %   N is the number of grid points
  %
  % Optional Inputs:
  %   alpha is the oversampling factor > 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Output:
  %   F the estimates of the Fourier coefficients along the trajectory
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on EE369C notes by John Pauly and Beatty et. al., IEEE TMI, 2005

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

  % % Density pre-compensation weights
  % weights = precomp_1D( traj, N, 'alpha', alpha, 'W', W, 'nC', nC );
  
  paddedData = iGridT_1D( F, traj, nGrid, ...
    'alpha', trueAlpha, 'W', W, 'nC', nC );

  % Crop out just the region of interest
  out = 1/N * cropData( paddedData, N );

end
