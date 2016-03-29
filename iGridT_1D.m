
function out = iGridT_1D( k, traj, N, varargin )
  % out = iGridT_1D( k, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   k is a 1D array of M elements specifying the k-space data values
  %   traj is a M element array specifying the k-space trajectory.
  %     The units are normalized to [-0.5,0.5).
  %   N is the number of elements of the output
  %
  % Optional Inputs:
  %   alpha is a float parameter specifying the oversampling factor
  %   W is an integer specifying the kernel's width
  %   nC specifies the number of samples in the kernel
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on EE369C notes written by John Pauly

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  nGrid = ceil( N * alpha );

  % Make the Kaiser Bessel convolution kernel
  G = nGrid;
  [kC,C,c1D,kw] = makeKbKernel( G, N, alpha, W, nC );

  gridKs = size2fftCoordinates( nGrid );

  nTraj = numel(traj);
  kDistThresh = 0.5*kw;
  kGridded = zeros( nGrid, 1 );
  LTraj = traj - 1;
  UTraj = traj + 1;
  for i=1:nTraj
    if mod(i,1000)==0
      fprintf('Working on %i of %i\n', i, nTraj );
    end
    kDists = abs( traj(i) - gridKs );
    shortDists = kDists( kDists < kDistThresh );
    CVals = interp1( kC, C, shortDists, 'linear', 0 );
    kGridded( kDists < kDistThresh ) = ...
     kGridded( kDists < kDistThresh ) + k(i) .* CVals;

    % LTraj and UTraj are used to accomplish circular convolution
    LkDists = abs( LTraj(i) - gridKs );
    LShortDists = LkDists( LkDists < kDistThresh );
    LCVals = interp1( kC, C, LShortDists, 'linear', 0 );
    kGridded( LkDists < kDistThresh ) = ...
      kGridded( LkDists < kDistThresh ) + k(i) .* LCVals;

    UkDists = abs( UTraj(i) - gridKs );
    UShortDists = UkDists( UkDists < kDistThresh );
    UCVals = interp1( kC, C, UShortDists, 'linear', 0 );
    kGridded( UkDists < kDistThresh ) = ...
      kGridded( UkDists < kDistThresh ) + k(i) .* UCVals;
  end

  data = alpha/nGrid * fftshift( ifft( ifftshift(kGridded) ) );
    % Divide by nGrid to account for ifft of convolution

  % Extract out the center region
  extracted = cropData( data, N );

  % Perform deapodization
  out = extracted ./ transpose(c1D);
end

