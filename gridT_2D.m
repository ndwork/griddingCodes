
function out = gridT_2D( in, traj, N, weights, varargin )
  % out = gridT_2D( in, traj, N, weights, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI reconstruction with Gridding
  %
  % Inputs:
  %   in is a 1D array representing the Fourier values
  %   traj is a Mx2 element array specifying the k-space trajectory.
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
  %   out is the uniformly spaced data in the space domain
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
  trueAlpha = max( nGrid ./ N );

  padded = padData(in, nGrid);
  tmp = iGrid_2D( padded, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
  out = tmp .* weights;
end

