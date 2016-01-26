
function out = iGridT_2D( k, traj, N, varargin )
  % out = iGridT_2D( k, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC,
  %   'verbose', verbose ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   k is a 1D array of M elements specifying the k-space data values
  %   traj is a Mx2 array specifying the k-space trajectory.
  %     The first/second column is kx/ky
  %     The units are normalized to [-0.5,0.5).
  %   N is the size of the output image
  %     If N is a scalar, then the final image is assumed to be square
  %     If N is a 2 element array, then N = [Ny Nx]
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
    Ny=N;  Nx=N; 
  else
    Ny=N(1); Nx=N(2);
  end

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

  nKy = ceil( Ny * alpha );
  nKx = ceil( Nx * alpha );

  % Make the convolution kernel
  G = max( nKy, nKx );
  kw = W/G;
  beta = pi * sqrt( W*W/(alpha*alpha) * (alpha-0.5)^2 - 0.8 );
  kC = linspace(0, 0.5*kw, nC);
  C = 1/kw * besseli( 0, beta * sqrt( 1 - ( 2*kC/kw ).^2 ) );

  gridKs = size2fftCoordinates([ nKy nKx ]);
  gridKy=gridKs{1};  gridKx=gridKs{2};
  [gridKx,gridKy] = meshgrid(gridKx,gridKy);

  nTraj = size(traj,1);
  kDistThresh = 0.5*kw;
  kGridded = zeros( nKy, nKx );
  UN_Traj = traj + [  ones(nTraj,1) zeros(nTraj,1) ];
  LN_Traj = traj + [ -ones(nTraj,1) zeros(nTraj,1) ];
  NU_Traj = traj + [  zeros(nTraj,1)  ones(nTraj,1) ];
  NL_Traj = traj + [  zeros(nTraj,1) -ones(nTraj,1) ];
  for trajIndx=1:nTraj
    % Note: this is convolution, but need to do circular convolution
    if verbose==true && mod(trajIndx,1000)==0
      fprintf('Working on %i of %i\n', trajIndx, nTraj );
    end
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThresh & distsKx < kDistThresh );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );

    CValsY = interp1( kC, C, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kC, C, shortDistsKx, 'linear', 0 );

    kGridded( shortDistIndxs ) = kGridded( shortDistIndxs ) + ...
      k(trajIndx) .* ( CValsY .* CValsX );

    % Implement circular convolution
    UN_kyDists = abs( UN_Traj(trajIndx,1) - gridKy );    % Upper -  Normal
    UN_kxDists = abs( UN_Traj(trajIndx,2) - gridKx );
    UN_shortDistIndxs = find( UN_kyDists < kDistThresh & UN_kxDists < kDistThresh );
    if numel( UN_shortDistIndxs ) > 0
      UN_shortDistsKy = UN_kyDists( UN_shortDistIndxs );
      UN_shortDistsKx = UN_kxDists( UN_shortDistIndxs );
      UN_CValsY = interp1( kC, C, UN_shortDistsKy, 'linear', 0 );
      UN_CValsX = interp1( kC, C, UN_shortDistsKx, 'linear', 0 );
      kGridded( UN_shortDistIndxs ) = kGridded( UN_shortDistIndxs ) + ...
        k(trajIndx) .* (UN_CValsY .* UN_CValsX);
    end

    LN_kyDists = abs( LN_Traj(trajIndx,1) - gridKy );    % Lower -  Normal
    LN_kxDists = abs( LN_Traj(trajIndx,2) - gridKx );
    LN_shortDistIndxs = find( LN_kyDists < kDistThresh & LN_kxDists < kDistThresh );
    if numel( LN_shortDistIndxs ) > 0
      LN_shortDistsKy = LN_kyDists( LN_shortDistIndxs );
      LN_shortDistsKx = LN_kxDists( LN_shortDistIndxs );
      LN_CValsY = interp1( kC, C, LN_shortDistsKy, 'linear', 0 );
      LN_CValsX = interp1( kC, C, LN_shortDistsKx, 'linear', 0 );
      kGridded( LN_shortDistIndxs ) = kGridded( LN_shortDistIndxs ) + ...
        k(trajIndx) .* (LN_CValsY .* LN_CValsX);
    end

    NU_kyDists = abs( NU_Traj(trajIndx,1) - gridKy );    % Normal - Upper
    NU_kxDists = abs( NU_Traj(trajIndx,2) - gridKx );
    NU_shortDistIndxs = find( NU_kyDists < kDistThresh & NU_kxDists < kDistThresh );
    if numel( NU_shortDistIndxs ) > 0
      NU_shortDistsKy = NU_kyDists( NU_shortDistIndxs );
      NU_shortDistsKx = NU_kxDists( NU_shortDistIndxs );
      NU_CValsY = interp1( kC, C, NU_shortDistsKy, 'linear', 0 );
      NU_CValsX = interp1( kC, C, NU_shortDistsKx, 'linear', 0 );
      kGridded( NU_shortDistIndxs ) = kGridded( NU_shortDistIndxs ) + ...
        k(trajIndx) .* (NU_CValsY .* NU_CValsX);
    end

    NL_kyDists = abs( NL_Traj(trajIndx,1) - gridKy );    % Normal - Lower
    NL_kxDists = abs( NL_Traj(trajIndx,2) - gridKx );
    NL_shortDistIndxs = find( ...
      NL_kyDists < kDistThresh & NL_kxDists < kDistThresh );
    if numel( NL_shortDistIndxs ) > 0
      NL_shortDistsKy = NL_kyDists( NL_shortDistIndxs );
      NL_shortDistsKx = NL_kxDists( NL_shortDistIndxs );
      NL_CValsY = interp1( kC, C, NL_shortDistsKy, 'linear', 0 );
      NL_CValsX = interp1( kC, C, NL_shortDistsKx, 'linear', 0 );
      kGridded( NL_shortDistIndxs ) = kGridded( NL_shortDistIndxs ) + ...
        k(trajIndx) .* (NL_CValsY .* NL_CValsX);
    end
  end

  img = alpha/nKy * alpha/nKx * fftshift( ifft2( ifftshift(kGridded) ) );

  % Extract out the center region
  extracted = cropImg( img, N );

  % Perform deapodization
  [y,x] = size2imgCoordinates( size( extracted ) );
  tmpX = sqrt( (pi*kw*x).^2 - beta*beta );
  cImgX = sinc( tmpX / pi );
  tmpY = sqrt( (pi*kw*y).^2 - beta*beta );
  cImgY = sinc( tmpY / pi );
  [cImgX,cImgY] = meshgrid( cImgX, cImgY );
  cImg = cImgX .* cImgY;
  out = extracted ./ cImg;
end

