
function out = iGridT_1D( F, traj, N, varargin )
  % out = iGridT_1D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
  %     along the trajectory.
  %   traj is a M element array specifying the k-space trajectory.
  %     The units are normalized to [-0.5,0.5).
  %   N is the number of elements of the output
  %
  % Optional Inputs:
  %   alpha is a float parameter specifying the oversampling factor
  %   W is an integer specifying the kernel's width
  %   nC specifies the number of samples in the kernel
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on EE369C notes written by John Pauly

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultVerbose = false;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.addParamValue( 'verbose', defaultVerbose );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;


  % Make the Kaiser Bessel convolution kernel
  G = N;
  [kC,C,c1D,kw] = makeKbKernel( G, N, alpha, W, nC );

  % Perform a circular convolution
  fftGridded = applyC_1D( F, traj, N, kw, kC, C );

  data = fftshift( ifft( ifftshift(fftGridded) ) );

  % Perform deapodization
  out = data ./ transpose(c1D);
end

