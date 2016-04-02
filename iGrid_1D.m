
function F = iGrid_1D( data, traj, varargin )
  % F = iGrid_1D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI encoding with Inverse Gridding
  %
  % Inputs
  %   data is a 1D array to be encoded
  %   traj is a M element array specifying the k-space trajectory.
  %     The units are normalized to [-0.5,0.5).
  %
  % Optional Inputs:
  %   alpha is the oversampling factor > 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Output:
  %   F the estimates of the Fourier coefficients along the trajectory
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on EE369C notes by John Pauly and Beatty et. al., IEEE TMI, 2005

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

  nData = numel(data);

  % Make the Kaiser Bessel convolution kernel
  nGrid = nData;
  G = nGrid;
  [kC,C,c1D,kw] = makeKbKernel( G, nGrid, alpha, W, nC );

  % Pre-emphasize the image
  preEmphasized = data ./ transpose(c1D);

  fftData = 1/nGrid * fftshift( fft( ifftshift(preEmphasized) ) );
    % divide my N to account for convolution

  % Perform a circular convolution
  gridKs = size2fftCoordinates( nGrid );
  nTraj = numel(traj);
  F = zeros( nTraj, 1 );
  kDistThresh = 0.5*kw;
  for trajIndx = 1:nTraj
    kDists = abs( traj(trajIndx) - gridKs );
    shortDistIndxs = find( kDists < kDistThresh );
    shortDists = kDists( shortDistIndxs );
    CVals = interp1( kC, C, shortDists, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    F(trajIndx) = F(trajIndx) + sum( kVals .* CVals );
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
      NewKVals = fftData( NewShortDistIndxs );
      F(trajIndx) = F(trajIndx) + sum( NewKVals .* NewCVals );
    end
  end

end

