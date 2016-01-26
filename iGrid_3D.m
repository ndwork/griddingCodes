
function k = applyE_3D( img, traj, varargin )
  % MRI encoding with Inverse Gridding
  % img is a 3D array specifying the volume to be encoded
  % traj is a Mx3 array specifying the k-space trajectory.
  %   The first/second/third column is kx/ky/kz
  %   The units are normalized to [-0.5,0.5).
  % k is a 1D array of M elements specifying the k-space data values

  defaultAlpha = 1.25;
  defaultW = 6;
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

  [Ny,Nx,Nz] = size( img );
  nKwy = ceil( Ny * alpha );
  nKwx = ceil( Nx * alpha );
  nKwz = ceil( Nz * alpha );

  % Make the convolution kernel
  G = max([ nKwy, nKwx, nKwz ]);
  beta = pi * sqrt( W*W/(alpha*alpha) * (alpha-0.5)^2 - 0.8 );
  ck = linspace(0, W/(2*G), nC);
  c = G/W * besseli( 0, beta * sqrt( 1 - ( 2*G*ck/W ).^2 ) );

  kwy = linspace( -0.5, 0.5, nKwy );
  kwx = linspace( -0.5, 0.5, nKwx );
  kwz = linspace( -0.5, 0.5, nKwz );
  [kwy,kwx,kwz] = meshgrid( kwy, kwx, kwz );

  kr = sqrt( kwy.*kwy + kwx.*kwx + kwz.*kwz );
  c = interp1( ck, c, kr, 'linear', 0 );
  ifftC = fftshift( ifftn( ifftshift(c) ) );

  padded = zeros( nKwy, nKwx, nKwz );
  minY = ceil( nKwy/2 - Ny/2 + 1 );
  minX = ceil( nKwx/2 - Nx/2 + 1 );
  minZ = ceil( nKwz/2 - Nz/2 + 1 );
  padded( minY:minY+Ny-1, minX:minX+Nx-1, minZ:minZ+Nz-1 ) = img;
  preEmphasized = padded ./ ifftC;
  fftPadded = fftshift( fftn( ifftshift(preEmphasized) ) );

  nK = size( traj, 1 );
  k = zeros( nK, 1 );
  for kIndx = 1:nK
    distsY = ( traj(kIndx,1) - kwy );  distsYSq = distsY .* distsY;
    distsX = ( traj(kIndx,2) - kwx );  distsXSq = distsX .* distsX;
    distsZ = ( traj(kIndx,3) - kwz );  distsZSq = distsZ .* distsZ;

    dists = sqrt( distsYSq + distsXSq + distsZSq );
    shortDists = dists( dists < kw/2 );
    cVals = interp1( ck, c, shortDists, 'linear', 0 );
    kVals = fftPadded( dists < kw/2 );
    k(kIndx) = sum( kVals .* cVals );
  end

end


