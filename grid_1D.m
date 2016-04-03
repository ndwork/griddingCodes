
function out = grid_1D( F, traj, N, varargin )

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

  nGrid = ceil( alpha * N );

  % Determine the density compensation weights
  
  % Make the Kaiser Bessel convolution kernel
  G = nGrid;
  [kC,C,c1D,kw] = makeKbKernel( G, nGrid, alpha, W, nC );

  % Perform a circular convolution
  gridKs = size2fftCoordinates( nGrid );
  nTraj = numel(traj);
  kDistThresh = 0.5*kw;
  fftData = zeros(nGrid,1);
  for trajIndx = 1:nTraj
    kDists = abs( traj(trajIndx) - gridKs );
    shortDistIndxs = find( kDists < kDistThresh );
    shortDists = kDists( shortDistIndxs );
    CVals = interp1( kC, C, shortDists, 'linear', 0 );
    fftData(shortDistIndxs) = fftData(shortDistIndxs) + ...
      F(trajIndx) * CVals;
  end

  for alt=[-1 1]
    NewTraj = traj + alt;
    if alt < 0
      NewTrajIndxs = find( NewTraj > -0.5-kw/2 );
    else
      NewTrajIndxs = find( NewTraj < 0.5+kw/2 );
    end
    NewTraj = NewTraj( NewTrajIndxs );

    for i=1:numel(NewTraj)
      trajIndx = NewTrajIndxs(i);
      NewkDists = abs( NewTraj(i) - gridKs );
      NewShortDistIndxs = find( NewkDists < kDistThresh );
      NewShortDists = NewkDists( NewShortDistIndxs );
      NewCVals = interp1( kC, C, NewShortDists, 'linear', 0 );
      fftData(NewShortDistIndxs) = fftData(NewShortDistIndxs) + ...
        F(trajIndx) * NewCVals;
    end
  end

  % Perform an Inverse DFT
  paddedData = alpha/nGrid * fftshift( ifft( ifftshift( fftData ) ) );

  % Crop out just the region of interest
  data = cropData( paddedData, N );
  croppedC1D = cropData( c1D, N );
  
  % Perform deapodization
  out = data ./ transpose(croppedC1D);
end
