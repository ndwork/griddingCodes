
function F = iGrid_1D( data, traj, varargin )
  % F = iGrid_1D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI encoding with Inverse Gridding
  %
  % Inputs
  %   data is a 1D array to be encoded
  %   traj is a M element array specifying the k-space trajectory.
  %     The units are normalized to [-0.5,0.5).
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

  N = numel(data);

  % Make the Kaiser Bessel convolution kernel
  G = N;
  [kC,C,c1D,kw] = makeKbKernel( G, N, alpha, W, nC );

  % Pre-emphasize the image
  preEmphasized = data ./ ( N * transpose(c1D) );

  % Perform DFT
  fftData = fftshift( fft( ifftshift(preEmphasized) ) );

  % Perform a circular convolution;
  F = applyCT_1D( fftData, traj, N, kw, kC, C );

end

