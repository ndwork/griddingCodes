
function fftTraj = iGrid_1D( data, traj, varargin )
  % k = iGrid_1D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % MRI encoding with Inverse Gridding
  %
  % Inputs
  %   data is a 1D array to be encoded
  %   traj is a M element array specifying the k-space trajectory.
  %   k is a 1D array of M elements specifying the k-space data values
  %     along the sepcified trajectory.
  %     The units are normalized to [-0.5,0.5).
  %
  % Optional Inputs:
  %   alpha is the oversampling factor > 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
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
  nGrid = ceil( nData * alpha );
  if mod(nData,2)==0
    minY = floor( nGrid/2 - nData/2 + 1 );
  else
    minY = ceil( nGrid/2 - nData/2 + 1 );
  end
  padded = zeros( nGrid, 1 );
  padded( minY : minY+nData-1 ) = data;

  % Make the Kaiser Bessel convolution kernel
  G = nGrid;
  [kC,C,c1D,kw] = makeKbKernel( G, nGrid, alpha, W, nC );

  % Pre-emphasize the image
  preEmphasized = padded ./ c1D';

  fftData = 1/nGrid * fftshift( fft( ifftshift(preEmphasized) ) );
    % Divide by nGrid to account for fft of convolution

  gridKs = size2fftCoordinates( nGrid );
  nTraj = numel(traj);
  fftTraj = zeros( nTraj, 1 );
  kDistThresh = 0.5*kw;
  LTraj = traj - 1;
  UTraj = traj + 1;
  for trajIndx = 1:nTraj
    kDists = abs( traj(trajIndx) - gridKs );
    shortDistIndxs = find( kDists < kDistThresh );
    shortDists = kDists( shortDistIndxs );
    CVals = interp1( kC, C, shortDists, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    fftTraj(trajIndx) = fftTraj(trajIndx) + sum( kVals .* CVals );

    % LTraj and UTraj are used to accomplish circular convolution
    LkDists = abs( LTraj(trajIndx) - gridKs );
    LShortDistIndxs = find( LkDists < kDistThresh );
    LShortDists = LkDists( LShortDistIndxs );
    LCVals = interp1( kC, C, LShortDists, 'linear', 0 );
    LKVals = fftData( LShortDistIndxs );
    fftTraj(trajIndx) = fftTraj(trajIndx) + sum( LKVals .* LCVals );

    UkDists = abs( UTraj(trajIndx) - gridKs );
    UShortDistIndxs = find( UkDists < kDistThresh );
    UShortDists = UkDists( UShortDistIndxs );
    UCVals = interp1( kC, C, UShortDists, 'linear', 0 );
    UKVals = fftData( UShortDistIndxs );
    fftTraj(trajIndx) = fftTraj(trajIndx) + sum( UKVals .* UCVals );
  end

end

