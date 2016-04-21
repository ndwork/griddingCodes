
function out = grid_1D( F, traj, N, weights, varargin )
  % F = grid_1D( F, traj, N, weights, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI reconstruction with Gridding
  %
  % Inputs:
  %   F is a 1D array representing the Fourier values
  %   traj is a M element array specifying the k-space trajectory.
  %     The units are normalized to [-0.5,0.5).
  %   N is the number of grid points
  %   weights is a 1D array; it is the pre-density compensation weights and
  %     can be generated using makePrecompWeights_1D.  Alternatively, they
  %     can be determined analytically for some sequences.
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

  weightedF = F .* weights;

  paddedData = iGridT_1D( weightedF, traj, nGrid, ...
    'alpha', trueAlpha, 'W', W, 'nC', nC );

  out = cropData( paddedData, N );
end
