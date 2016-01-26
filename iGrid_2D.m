
function kTraj = iGrid_2D( data, traj, varargin )
  % k = applyE_2D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI encoding with Inverse Gridding
  %
  % Inputs
  %   data is a 2D array specifying the volume to be encoded
  %   traj is a Mx2 array specifying the k-space trajectory.
  %     The first/second column is kx/ky
  %   k is a 1D array of M elements specifying the k-space data values
  %     The units are normalized to [-0.5,0.5).
  %
  % Optional Inputs:
  %   alpha is the oversampling factor >= 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on Beatty et. al., IEEE TMI, 2005

  defaultAlpha = 1.25;
  defaultW = 6;
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
  nGridY = ceil( Ny * alpha );
  nGridX = ceil( Nx * alpha );
  padded = zeros( nGridY, nGridX );
  minY = ceil( nGridY/2 - Ny/2 + 1 );
  minX = ceil( nGridX/2 - Nx/2 + 1 );
  padded( minY:minY+Ny-1, minX:minX+Nx-1 ) = data;

  % Pre-emphasize the image  
  beta = pi * sqrt( W*W/(alpha*alpha) * (alpha-0.5)^2 - 0.8 );
  G = max( nGridY, nGridX );
  kw = W/G;
  [y,x] = size2imgCoordinates( [nGridY nGridX] );
  tmpY = sqrt( (pi*kw*y).^2 - beta*beta );
  cImgY = kw * sinc( tmpY / pi );
  tmpX = sqrt( (pi*kw*x).^2 - beta*beta );
  cImgX = kw * sinc( tmpX / pi );
  [cImgX,cImgY] = meshgrid( cImgX, cImgY );
  cImg = cImgX .* cImgY;
  preEmphasized = padded ./ cImg;
  fftData = 1/(nGridY*nGridX) * fftshift( fft2( ifftshift(preEmphasized) ) );

  % Make the convolution kernel
  kC = linspace(0, 0.5*kw, nC)';
  C = besseli( 0, beta * sqrt( 1 - ( 2*kC/kw ).^2 ) );

  gridKs = size2fftCoordinates( [nGridY nGridX] );
  gridKys=gridKs{1};  gridKxs=gridKs{2};
  [gridKxs,gridKys] = meshgrid( gridKxs, gridKys );

  nTraj = size( traj, 1 );
  kDistThresh = 0.5*kw;
  kTraj = zeros( nTraj, 1 );
  UN_Traj = traj + [  ones(nTraj,1) zeros(nTraj,1) ];
  LN_Traj = traj + [ -ones(nTraj,1) zeros(nTraj,1) ];
  NU_Traj = traj + [  zeros(nTraj,1)  ones(nTraj,1) ];
  NL_Traj = traj + [  zeros(nTraj,1) -ones(nTraj,1) ];
  for trajIndx = 1:nTraj
    distsKy = ( traj(trajIndx,1) - gridKys );
    distsKx = ( traj(trajIndx,2) - gridKxs );
    shortDistIndxs = find( distsKy < kDistThresh & distsKx < kDistThresh );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kC, C, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kC, C, shortDistsKx, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    kTraj( trajIndx ) = kTraj( trajIndx ) + sum( kVals .* CValsY .* CValsX );

    % Implement circular convolution
    UN_distsKy = ( UN_Traj(trajIndx,1) - gridKys );
    UN_distsKx = ( UN_Traj(trajIndx,2) - gridKxs );
    UN_shortDistIndxs = find( UN_distsKy < kDistThresh & UN_distsKx < kDistThresh );
    UN_shortDistsKy = UN_distsKy( UN_shortDistIndxs );
    UN_shortDistsKx = UN_distsKx( UN_shortDistIndxs );
    UN_CValsY = interp1( kC, C, UN_shortDistsKy, 'linear', 0 );
    UN_CValsX = interp1( kC, C, UN_shortDistsKx, 'linear', 0 );
    UN_kVals = fftData( UN_shortDistIndxs );
    kTraj( trajIndx ) = kTraj( trajIndx ) + sum( UN_kVals .* UN_CValsY .* UN_CValsX );

    LN_distsKy = ( LN_Traj(trajIndx,1) - gridKys );
    LN_distsKx = ( LN_Traj(trajIndx,2) - gridKxs );
    LN_shortDistIndxs = find( LN_distsKy < kDistThresh & LN_distsKx < kDistThresh );
    LN_shortDistsKy = LN_distsKy( LN_shortDistIndxs );
    LN_shortDistsKx = LN_distsKx( LN_shortDistIndxs );
    LN_CValsY = interp1( kC, C, LN_shortDistsKy, 'linear', 0 );
    LN_CValsX = interp1( kC, C, LN_shortDistsKx, 'linear', 0 );
    LN_kVals = fftData( LN_shortDistIndxs );
    kTraj( trajIndx ) = kTraj( trajIndx ) + sum( LN_kVals .* LN_CValsY .* LN_CValsX );

    NU_distsKy = ( NU_Traj(trajIndx,1) - gridKys );
    NU_distsKx = ( NU_Traj(trajIndx,2) - gridKxs );
    NU_shortDistIndxs = find( NU_distsKy < kDistThresh & NU_distsKx < kDistThresh );
    NU_shortDistsKy = distsKy( NU_shortDistIndxs );
    NU_shortDistsKx = distsKx( NU_shortDistIndxs );
    NU_CValsY = interp1( kC, C, NU_shortDistsKy, 'linear', 0 );
    NU_CValsX = interp1( kC, C, NU_shortDistsKx, 'linear', 0 );
    NU_kVals = fftData( NU_shortDistIndxs );
    kTraj( trajIndx ) = kTraj( trajIndx ) + sum( NU_kVals .* NU_CValsY .* NU_CValsX );

    NL_distsKy = ( NL_Traj(trajIndx,1) - gridKys );
    NL_distsKx = ( NL_Traj(trajIndx,2) - gridKxs );
    NL_shortDistIndxs = find( NL_distsKy < kDistThresh & NL_distsKx < kDistThresh );
    NL_shortDistsKy = distsKy( NL_shortDistIndxs );
    NL_shortDistsKx = distsKx( NL_shortDistIndxs );
    NL_CValsY = interp1( kC, C, NL_shortDistsKy, 'linear', 0 );
    NL_CValsX = interp1( kC, C, NL_shortDistsKx, 'linear', 0 );
    NL_kVals = fftData( NL_shortDistIndxs );
    kTraj( trajIndx ) = kTraj( trajIndx ) + sum( NL_kVals .* NL_CValsY .* NL_CValsX );
  end

end
