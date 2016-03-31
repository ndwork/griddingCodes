
function out = iGridT_2D( F, traj, N, varargin )
  % out = iGridT_2D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC,
  %   'verbose', verbose ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
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

  %nGridY = ceil( Ny * alpha );
  %nGridX = ceil( Nx * alpha );
  nGridY = Ny;
  nGridX = Nx;

  %% Make the convolution kernel
  Gy = nGridY;
  [kCy,Cy,cImgY,kwy] = makeKbKernel( Gy, Ny, alpha, W, nC );
  Gx = nGridX;
  [kCx,Cx,cImgX,kwx] = makeKbKernel( Gx, Nx, alpha, W, nC );

  gridKs = size2fftCoordinates([ nGridY nGridX ]);
  gridKy=gridKs{1};  gridKx=gridKs{2};
  [gridKx,gridKy] = meshgrid(gridKx,gridKy);

  % Perform Adjoint of Circular Convolution
  nTraj = size(traj,1);
  kDistThreshY = 0.5*kwy;
  kDistThreshX = 0.5*kwx;
  fftGridded = zeros( nGridY, nGridX );
  UN_traj = traj + [  ones(nTraj,1) zeros(nTraj,1) ];
  LN_traj = traj + [ -ones(nTraj,1) zeros(nTraj,1) ];
  NU_traj = traj + [  zeros(nTraj,1)  ones(nTraj,1) ];
  NL_traj = traj + [  zeros(nTraj,1) -ones(nTraj,1) ];
  for trajIndx=1:nTraj
    if verbose==true && mod(trajIndx,1000)==0
      fprintf('Working on %i of %i\n', trajIndx, nTraj );
    end

    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThreshY & distsKx < kDistThreshX );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    fftGridded(shortDistIndxs) = fftGridded(shortDistIndxs) + ...
      F(trajIndx) * ( CValsY .* CValsX );

    UN_distsKy = abs( UN_traj(trajIndx,1) - gridKy );
    UN_distsKx = abs( UN_traj(trajIndx,2) - gridKx );
    UN_shortDistIndxs = find( UN_distsKy < kDistThreshY & UN_distsKx < kDistThreshX );
    UN_shortDistsKy = UN_distsKy( UN_shortDistIndxs );
    UN_shortDistsKx = UN_distsKx( UN_shortDistIndxs );
    UN_CValsY = interp1( kCy, Cy, UN_shortDistsKy, 'linear', 0 );
    UN_CValsX = interp1( kCx, Cx, UN_shortDistsKx, 'linear', 0 );
    fftGridded(UN_shortDistIndxs) = fftGridded(UN_shortDistIndxs) + ...
      F(trajIndx) * ( UN_CValsY .* UN_CValsX );

    LN_distsKy = abs( LN_traj(trajIndx,1) - gridKy );
    LN_distsKx = abs( LN_traj(trajIndx,2) - gridKx );
    LN_shortDistIndxs = find( LN_distsKy < kDistThreshY & LN_distsKx < kDistThreshX );
    LN_shortDistsKy = LN_distsKy( LN_shortDistIndxs );
    LN_shortDistsKx = LN_distsKx( LN_shortDistIndxs );
    LN_CValsY = interp1( kCy, Cy, LN_shortDistsKy, 'linear', 0 );
    LN_CValsX = interp1( kCx, Cx, LN_shortDistsKx, 'linear', 0 );
    fftGridded(LN_shortDistIndxs) = fftGridded(LN_shortDistIndxs) + ...
      F(trajIndx) * ( LN_CValsY .* LN_CValsX );

    NU_distsKy = abs( NU_traj(trajIndx,1) - gridKy );
    NU_distsKx = abs( NU_traj(trajIndx,2) - gridKx );
    NU_shortDistIndxs = find( NU_distsKy < kDistThreshY & NU_distsKx < kDistThreshX );
    NU_shortDistsKy = NU_distsKy( NU_shortDistIndxs );
    NU_shortDistsKx = NU_distsKx( NU_shortDistIndxs );
    NU_CValsY = interp1( kCy, Cy, NU_shortDistsKy, 'linear', 0 );
    NU_CValsX = interp1( kCx, Cx, NU_shortDistsKx, 'linear', 0 );
    fftGridded(NU_shortDistIndxs) = fftGridded(NU_shortDistIndxs) + ...
      F(trajIndx) * ( NU_CValsY .* NU_CValsX );

    NL_distsKy = abs( NL_traj(trajIndx,1) - gridKy );
    NL_distsKx = abs( NL_traj(trajIndx,2) - gridKx );
    NL_shortDistIndxs = find( NL_distsKy < kDistThreshY & NL_distsKx < kDistThreshX );
    NL_shortDistsKy = NL_distsKy( NL_shortDistIndxs );
    NL_shortDistsKx = NL_distsKx( NL_shortDistIndxs );
    NL_CValsY = interp1( kCy, Cy, NL_shortDistsKy, 'linear', 0 );
    NL_CValsX = interp1( kCx, Cx, NL_shortDistsKx, 'linear', 0 );
    fftGridded(NL_shortDistIndxs) = fftGridded(NL_shortDistIndxs) + ...
      F(trajIndx) * ( NL_CValsY .* NL_CValsX );
  end

  % Perform an inverse fft
  data = fftshift( ifft2( ifftshift(fftGridded) ) );

  % Extract out the center region
  %extracted = cropImg( img, N );

  % Perform deapodization
  cImg = transpose(cImgY) * cImgX;
  %out = extracted ./ cImg;
  out = data ./ cImg;
end

