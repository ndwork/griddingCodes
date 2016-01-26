
function out = iGridT_3D( k, traj, N, varargin )
  % out = applyET( k, traj, N, [ 'alpha', alpha, 'w', w, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % (often called inverse gridding).  This function applies the transpose 
  % of inverse gridding to the input data.
  % k is a 1D array of M elements specifying the k-space data values
  % traj is a Mx3 array specifying the k-space trajectory.
  %   The first/second/third column is kx/ky/kz
  %   The units are normalized to [-0.5,0.5).
  % N is the size of the output image
  %   If N is a scalar, then the final image is assumed to be square
  %   If N is a 3 element array, then N = [Nx Ny Nz]
  % alpha is an optional float parameter specifying the oversampling factor
  % w is an integer specifying the kernel's width
  % nC specifies the number of samples in the kernel

  if numel(N)==1
    Nx=N; Ny=N; Nz=N;
  else
    Ny=N(1); Nx=N(2); Nz=N(3);
  end

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'w', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  w = p.Results.w;
  nC = p.Results.nC;

  % Make the convolution kernel
  G = max([ nKwy, nKwx, nKwz ]);
  beta = pi * sqrt( W*W/(alpha*alpha) * (alpha-0.5)^2 - 0.8 );
  ck = linspace(0, W/(2*G), nC);
  c = G/W * besseli( 0, beta * sqrt( 1 - ( 2*G*ck/W ).^2 ) );

  nKwy = ceil( Ny * alpha );  kwy = linspace( -0.5, 0.5, nKwy );
  nKwx = ceil( Nx * alpha );  kwx = linspace( -0.5, 0.5, nKwx );
  nKwz = ceil( Nz * alpha );  kwz = linspace( -0.5, 0.5, nKwz );

  kGridded = zeros( nKwx, nKwy, nKwz );
  for kwyIndx=1:nKwy
    distsYsq = ( traj(:,1) - kwy(kwyIndx) ).^2;

    for kwxIndx=1:nKwx
      distsXsq = ( traj(:,2) - kwx(kwxIndx) ).^2;

      for kwzIndx=1:nKwz
        distsZsq = ( traj(:,3) - kwz(kwzIndx) ).^2;

        dists = sqrt( distsXsq + distsYsq + distsZsq );
        shortDists = dists( dists < kw/2 );
        kVals = k( dists < kw/2 );
        cVals = interp1( ck, c, shortDists, 'linear', 0 );
        kGridValues = kVals .* cVals;

        kGridded( kwyIndx, kwxIndx, kwzIndx ) = sum( kGridValues );
      end
    end
  end
  img = fftshift( ifftn( ifftshift(kGridded) ) );

  % Perform deapodization
  y = linspace( -0.5*alpha, 0.5*alpha, ceil(Nx*alpha) );
  x = linspace( -0.5*alpha, 0.5*alpha, ceil(Ny*alpha) );
  z = linspace( -0.5*alpha, 0.5*alpha, ceil(Nz*alpha) );
  [x,y,z] = meshgrid( x, y, z );
  r = sqrt( x.*x + y.*y + z.*z );
  tmp = sqrt( (pi*kw*r).^2 - beta^2 ) / pi;
  cImg = sinc( tmp );
  deapImg = img ./ cImg;

  % Extract out the center region
  sDeapImg = size( deapImg );
  minY = ceil( sDeapImg(1)/2 - Ny/2 + 1 );
  minX = ceil( sDeapImg(2)/2 - Nx/2 + 1 );
  minZ = ceil( sDeapImg(3)/2 - Nz/2 + 1 );
  out = deapImg( minY:minY+Ny-1, minX:minX+Nx-1, minZ:minZ+Nz-1 );
end




