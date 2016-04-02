
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
  kws = [ kwy kwx ];
  nGridZ = Nz;
  Gz = nGridZ;
  [kCz,Cz,cImgZ,kwz] = makeKbKernel( Gz, nGridZ, alpha, W, nC );
  kws = [ kwy, kwx, kwz ];

  % Pre-emphasize the image
  cImgYX = transpose(cImgY) * cImgX;
  cImgZ_reshaped = reshape( cImgZ, [1 1 numel(cImgZ)] );
  cImg = bsxfun( @times, cImgZ_reshaped, cImgYX );
  preEmphasized = data ./ cImg;

  % Perform an fft
  fftData = 1/(nGridY*nGridX*nGridZ) * fftshift( fftn( ifftshift(preEmphasized) ) );

  % Perform a circular convolution
  gridKs = size2fftCoordinates( [nGridY nGridX nGridZ] );
  gridKy=gridKs{1};  gridKx=gridKs{2};  gridKz=gridKs{3};
  [gridKx,gridKy,gridKz] = meshgrid( gridKx, gridKy, gridKz );

  nTraj = size( traj, 1 );
  kDistThreshY = 0.5*kwy;
  kDistThreshX = 0.5*kwx;
  kDistThreshZ = 0.5*kwz;
  F = zeros( nTraj, 1 );
  for trajIndx = 1:nTraj
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    distsKz = abs( traj(trajIndx,3) - gridKz );
    shortDistIndxs = find( distsKy < kDistThreshY & ...
                           distsKx < kDistThreshX & ...
                           distsKz < kDistThreshZ );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    shortDistsKz = distsKz( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    CValsZ = interp1( kCz, Cz, shortDistsKz, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    F( trajIndx ) = F( trajIndx ) + sum( kVals .* ...
      CValsY .* CValsX .* CValsZ );
  end

  % Circular convolution
  onesCol = ones(nTraj,1);
  for dim=1:3
    alt = zeros( size(traj) );
    alt(:,dim) = onesCol;

    for altDir=[-1 1]
      NewTraj = traj + altDir*alt;
      if altDir < 0
        NewTrajIndxs = find( NewTraj(:,dim) > -0.5-kws(dim)/2 );
      else
        NewTrajIndxs = find( NewTraj(:,dim) < 0.5+kws(dim)/2 );
      end

      NewTraj = NewTraj( NewTrajIndxs, : );
      for i=1:numel(NewTrajIndxs)
        trajIndx = NewTrajIndxs(i);
        NewDistsKy = abs( NewTraj(i,1) - gridKy );
        NewDistsKx = abs( NewTraj(i,2) - gridKx );
        NewDistsKz = abs( NewTraj(i,3) - gridKz );
        NewShortDistIndxs = find( NewDistsKy < kDistThreshY & ...
                                  NewDistsKx < kDistThreshX & ...
                                  NewDistsKz < kDistThreshZ );
        NewShortDistsKy = NewDistsKy( NewShortDistIndxs );
        NewShortDistsKx = NewDistsKx( NewShortDistIndxs );
        NewShortDistsKz = NewDistsKz( NewShortDistIndxs );
        NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
        NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
        NewCValsZ = interp1( kCz, Cz, NewShortDistsKz, 'linear', 0 );
        NewKVals = fftData( NewShortDistIndxs );
        F(trajIndx) = F(trajIndx) + sum( NewKVals .* ...
          NewCValsY .* NewCValsX .* NewCValsZ );
      end
    end
  end

end

