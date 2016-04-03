
function out = iGridT_1D( F, traj, N, varargin )
  % out = iGridT_1D( k, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
  %     along the trajectory.
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

  %nGrid = ceil( N * alpha );
  nGrid = N;

  % Make the Kaiser Bessel convolution kernel
  G = nGrid;
  [kC,C,c1D,kw] = makeKbKernel( G, N, alpha, W, nC );

  gridKs = size2fftCoordinates( nGrid );

  nTraj = numel(traj);
  kDistThresh = 0.5*kw;
  fftGridded = zeros( nGrid, 1 );
  LTraj = traj - 1;
  UTraj = traj + 1;
  for trajIndx=1:nTraj
    if verbose==true && mod(trajIndx,1000)==0
      fprintf('Working on %i of %i\n', trajIndx, nTraj );
    end

    kDists = abs( traj(trajIndx) - gridKs );
    shortDistIndxs = find( kDists < kDistThresh );
    shortDists = kDists( shortDistIndxs );
    CVals = interp1( kC, C, shortDists, 'linear', 0 );
    fftGridded(shortDistIndxs) = fftGridded(shortDistIndxs) + ...
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
      fftGridded(NewShortDistIndxs) = fftGridded(NewShortDistIndxs) + ...
        F(trajIndx) * NewCVals;
    end
  end

  data = fftshift( ifft( ifftshift(fftGridded) ) );
    % Divide by nGrid to account for ifft of convolution

  % Perform deapodization
  out = data ./ transpose(c1D);
end

