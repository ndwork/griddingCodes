
function recon = grid_2D( F, kTraj, N, weights, varargin )
  % recon = grid_2D( F, kTraj, N, weights, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Image reconstruction with Gridding
  %
  % Inputs:
  %   F is a 1D array representing the Fourier values
  %   kTraj is a Mx2 element array specifying the k-space trajectory.
  %     The first/second column are the kx/ky locations.
  %     The units are normalized to [-0.5,0.5).
  %   N is a 2 element array [Ny Nx] representing the number of grid points
  %   weights is a 1D array; it is the pre-density compensation weights and
  %     can be generated using makePrecompWeights_2D.  Alternatively, they
  %     can be determined analytically for some sequences.
  %
  % Optional Inputs:
  %   alpha is the oversampling factor > 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Output:
  %   recon is the uniformly spaced data in the space domain
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on EE369C notes by John Pauly and Beatty et. al., IEEE TMI, 2005

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

  weightedF = F .* weights;

  padded = iGridT_2D( weightedF, kTraj, nGrid, ...
    'alpha', trueAlpha, 'W', W, 'nC', nC );

  recon = cropData( padded, N );
end

