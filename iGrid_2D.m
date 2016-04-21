
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

  % Make the convolution kernel
  nGridY = Ny;
  Gy = nGridY;
  [kCy,Cy,cImgY,kwy] = makeKbKernel( Gy, nGridY, alpha, W, nC );
  nGridX = Nx;
  Gx = nGridX;
  [kCx,Cx,cImgX,kwx] = makeKbKernel( Gx, nGridX, alpha, W, nC );
  kws = [ kwy kwx ];

  % Pre-emphasize the image
  cImg = transpose(cImgY) * cImgX;
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
  for trajIndx = 1:nTraj
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThreshY & distsKx < kDistThreshX );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    F( trajIndx ) = sum( kVals .* CValsY .* CValsX );
  end

  % Circular convolution
  onesCol = ones(nTraj,1);
  for dim=1:2
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
        NewShortDistIndxs = find( NewDistsKy < kDistThreshY & ...
                                  NewDistsKx < kDistThreshX );
        NewShortDistsKy = NewDistsKy( NewShortDistIndxs );
        NewShortDistsKx = NewDistsKx( NewShortDistIndxs );
        NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
        NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
        NewKVals = fftData( NewShortDistIndxs );
        F(trajIndx) = F(trajIndx) + sum( NewKVals .* NewCValsY .* NewCValsX );
      end
    end
  end
  
end

