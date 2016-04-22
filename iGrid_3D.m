
function F = iGrid_3D( data, traj, varargin )
  % F = iGrid_3D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI encoding with Inverse Gridding
  %
  % Inputs
  %   data is a 2D array specifying the volume to be encoded
  %   traj is a Mx3 array specifying the k-space trajectory.
  %     The first/second/third column is ky/kx/kz
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
  % Written by Nicholas Dwork (c) 2015
  % Based on Beatty et. al., IEEE TMI, 2005

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

  [Ny,Nx,Nz] = size( data );

  % Make the convolution kernel
  nGridY = Ny;
  Gy = nGridY;
  [kCy,Cy,cImgY,kwy] = makeKbKernel( Gy, nGridY, alpha, W, nC );
  nGridX = Nx;
  Gx = nGridX;
  [kCx,Cx,cImgX,kwx] = makeKbKernel( Gx, nGridX, alpha, W, nC );
  nGridZ = Nz;
  Gz = nGridZ;
  [kCz,Cz,cImgZ,kwz] = makeKbKernel( Gz, nGridZ, alpha, W, nC );
  kws = [ kwy, kwx, kwz ];

  % Pre-emphasize the image
  cImgYX = (Ny*cImgY) * transpose(Nx*cImgX);
  cImgZ_reshaped = Nz * reshape( cImgZ, [1 1 numel(cImgZ)] );
  cImg = bsxfun( @times, cImgZ_reshaped, cImgYX );
  preEmphasized = data ./ cImg;

  % Perform an fft
  fftData = fftshift( fftn( ifftshift(preEmphasized) ) );

  % Perform a circular convolution
  N = [Ny Nx Nz];
  F = applyCT_3D( fftData, traj, N, kws, ...
    kCy, kCx, kCz, Cy, Cx, Cz );

end

