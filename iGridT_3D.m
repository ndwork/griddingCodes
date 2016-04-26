
function out = iGridT_3D( F, traj, N, varargin )
  % out = iGridT_3D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
  %   traj is a Mx3 array specifying the k-space trajectory.
  %     The first/second/third column is kx/ky/kz
  %     The units are normalized to [-0.5,0.5).
  %   N is the size of the output image
  %     If N is a scalar, then the final image is assumed to be square
  %     If N is a 2 element array, then N = [Ny Nx Nz]
  %
  % Optional Inputs:
  %   alpha - a float parameter specifying the oversampling factor
  %   W - an integer specifying the kernel's width
  %   nC - specifies the number of samples in the kernel
  %   verbose - if set to true, outputs processing status
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on EE369C notes written by John Pauly

  if numel(N)==1
    Ny=N;  Nx=N;  Nz=N;
  else
    Ny=N(1); Nx=N(2); Nz=N(3);
  end

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  % Make the convolution kernel
  Gy = Ny;
  [kCy,Cy,cImgY,kwy] = makeKbKernel( Gy, Ny, alpha, W, nC );
  Gx = Nx;
  [kCx,Cx,cImgX,kwx] = makeKbKernel( Gx, Nx, alpha, W, nC );
  Gz = Nz;
  [kCz,Cz,cImgZ,kwz] = makeKbKernel( Gz, Nz, alpha, W, nC );
  kws = [ kwy kwx kwz ];

  % Perform Adjoint of Circular Convolution
  fftGridded = applyC_3D( F, traj, [Ny Nx Nz], kws, ...
    kCy, kCx, kCz, Cy, Cx, Cz );

  % Perform an inverse fft
  data = fftshift( ifftn( ifftshift(fftGridded) ) );

  % Perform deapodization
  cImgYX = cImgY * transpose(cImgX);
  cImgZ_reshaped = reshape( cImgZ, [1 1 numel(cImgZ)] );
  cImg = bsxfun( @times, cImgZ_reshaped, cImgYX );
  out = data ./ cImg;
end

