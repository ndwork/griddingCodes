
function F = iGrid_2D( data, traj, varargin )
  % k = applyE_2D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
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

  [Ny,Nx] = size( data );
  %nGridY = ceil( Ny * alpha );
  %nGridX = ceil( Nx * alpha );
  %padded = zeros( nGridY, nGridX );
  %minY = ceil( nGridY/2 - Ny/2 + 1 );
  %minX = ceil( nGridX/2 - Nx/2 + 1 );
  %padded( minY:minY+Ny-1, minX:minX+Nx-1 ) = data;

  % Make the convolution kernel
  nGridY = Ny;
  Gy = nGridY;
  [kCy,Cy,cImgY,kwy] = makeKbKernel( Gy, nGridY, alpha, W, nC );
  nGridX = Nx;
  Gx = nGridX;
  [kCx,Cx,cImgX,kwx] = makeKbKernel( Gx, nGridX, alpha, W, nC );

  % Pre-emphasize the image
  cImg = transpose(cImgY) * cImgX;
  %preEmphasized = padded ./ cImg;
  preEmphasized = data ./ cImg;

  % Perform an fft
  fftData = 1/(nGridY*nGridX) * fftshift( fft2( ifftshift(preEmphasized) ) );

  % Perform a circular convolution
  gridKs = size2fftCoordinates( [nGridY nGridX] );
  gridKy=gridKs{1};  gridKx=gridKs{2};
  [gridKx,gridKy] = meshgrid( gridKx, gridKy );

  nTraj = size( traj, 1 );
  kDistThreshY = 0.5*kwy;
  kDistThreshX = 0.5*kwx;
  F = zeros( nTraj, 1 );
  UN_Traj = traj + [  ones(nTraj,1) zeros(nTraj,1) ];
  LN_Traj = traj + [ -ones(nTraj,1) zeros(nTraj,1) ];
  NU_Traj = traj + [  zeros(nTraj,1)  ones(nTraj,1) ];
  NL_Traj = traj + [  zeros(nTraj,1) -ones(nTraj,1) ];
  for trajIndx = 1:nTraj
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThreshY & distsKx < kDistThreshX );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    F( trajIndx ) = F( trajIndx ) + sum( kVals .* CValsY .* CValsX );

    UN_distsKy = abs( UN_Traj(trajIndx,1) - gridKy );
    UN_distsKx = abs( UN_Traj(trajIndx,2) - gridKx );
    UN_shortDistIndxs = find( UN_distsKy < kDistThreshY & UN_distsKx < kDistThreshX );
    UN_shortDistsKy = UN_distsKy( UN_shortDistIndxs );
    UN_shortDistsKx = UN_distsKx( UN_shortDistIndxs );
    UN_CValsY = interp1( kCy, Cy, UN_shortDistsKy, 'linear', 0 );
    UN_CValsX = interp1( kCx, Cx, UN_shortDistsKx, 'linear', 0 );
    UN_kVals = fftData( UN_shortDistIndxs );
    F( trajIndx ) = F( trajIndx ) + sum( UN_kVals .* UN_CValsY .* UN_CValsX );

    LN_distsKy = abs( LN_Traj(trajIndx,1) - gridKy );
    LN_distsKx = abs( LN_Traj(trajIndx,2) - gridKx );
    LN_shortDistIndxs = find( LN_distsKy < kDistThreshY & LN_distsKx < kDistThreshX );
    LN_shortDistsKy = LN_distsKy( LN_shortDistIndxs );
    LN_shortDistsKx = LN_distsKx( LN_shortDistIndxs );
    LN_CValsY = interp1( kCy, Cy, LN_shortDistsKy, 'linear', 0 );
    LN_CValsX = interp1( kCx, Cx, LN_shortDistsKx, 'linear', 0 );
    LN_kVals = fftData( LN_shortDistIndxs );
    F( trajIndx ) = F( trajIndx ) + sum( LN_kVals .* LN_CValsY .* LN_CValsX );

    NU_distsKy = abs( NU_Traj(trajIndx,1) - gridKy );
    NU_distsKx = abs( NU_Traj(trajIndx,2) - gridKx );
    NU_shortDistIndxs = find( NU_distsKy < kDistThreshY & NU_distsKx < kDistThreshX );
    NU_shortDistsKy = NU_distsKy( NU_shortDistIndxs );
    NU_shortDistsKx = NU_distsKx( NU_shortDistIndxs );
    NU_CValsY = interp1( kCy, Cy, NU_shortDistsKy, 'linear', 0 );
    NU_CValsX = interp1( kCx, Cx, NU_shortDistsKx, 'linear', 0 );
    NU_kVals = fftData( NU_shortDistIndxs );
    F( trajIndx ) = F( trajIndx ) + sum( NU_kVals .* NU_CValsY .* NU_CValsX );

    NL_distsKy = abs( NL_Traj(trajIndx,1) - gridKy );
    NL_distsKx = abs( NL_Traj(trajIndx,2) - gridKx );
    NL_shortDistIndxs = find( NL_distsKy < kDistThreshY & NL_distsKx < kDistThreshX );
    NL_shortDistsKy = NL_distsKy( NL_shortDistIndxs );
    NL_shortDistsKx = NL_distsKx( NL_shortDistIndxs );
    NL_CValsY = interp1( kCy, Cy, NL_shortDistsKy, 'linear', 0 );
    NL_CValsX = interp1( kCx, Cx, NL_shortDistsKx, 'linear', 0 );
    NL_kVals = fftData( NL_shortDistIndxs );
    F( trajIndx ) = F( trajIndx ) + sum( NL_kVals .* NL_CValsY .* NL_CValsX );
  end

end

