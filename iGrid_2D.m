
function F = iGrid_2D( data, traj, varargin )
  % F = iGrid_2D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI encoding with Inverse Gridding
  %
  % Inputs
  %   data is a 2D array specifying the volume to be encoded
  %   traj is a Mx2 array specifying the k-space trajectory.
  %     The first/second column is ky/kx
  %     The units are normalized to [-0.5,0.5).
  %
  % Optional Inputs:
  %   alpha is the oversampling factor >= 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Output:
  %   F the estimates of the Fourier coefficients along the trajectory
  %
  % Written by Nicholas Dwork (c) 2016
  % Based on EE369C notes by John Pauly and Beatty et. al., IEEE TMI, 2005

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  [Ny,Nx] = size( data );

  % Make the Kaiser Bessel convolution kernel
  Gy = Ny;
  [kCy,Cy,cImgY] = makeKbKernel( Gy, Ny, alpha, W, nC );
  Gx = Nx;
  [kCx,Cx,cImgX] = makeKbKernel( Gx, Nx, alpha, W, nC );

  % Pre-emphasize the image
  denom = cImgY * transpose(cImgX);
  preEmphasized = data ./ denom;

  % Perform an fft
  fftData = fftshift( fft2( ifftshift(preEmphasized) ) ) / (Ny*Nx);

  % Perform a circular convolution
  N = [Ny Nx];
  F = applyCT_2D( fftData, traj, N, kCy, kCx, Cy, Cx );
end

