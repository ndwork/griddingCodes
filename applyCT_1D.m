
function out = applyCT_1D( fftData, traj, N, kC, C )
  % out = applyCT_1D( fftData, traj, N, kC, C )
  %
  % Written by Nicholas Dwork - Copyright 2016

  gridKs = size2fftCoordinates( N );

  kw = max(kC);
  kDistThresh = 0.5*kw;

  nTraj = numel( traj );
  out = zeros( nTraj, 1 );  
  for trajIndx = 1:nTraj
    kDists = abs( traj(trajIndx) - gridKs );
    shortDistIndxs = find( kDists < kDistThresh );
    shortDists = kDists( shortDistIndxs );
    CVals = interp1( kC, C, shortDists, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    out(trajIndx) = sum( kVals .* CVals );
  end

  for alt=[-1 1]
    NewTraj = traj + alt;
    if alt < 0
      NewTrajIndxs = find( NewTraj > -0.5-kw );
    else
      NewTrajIndxs = find( NewTraj < 0.5+kw );
    end
    NewTraj = NewTraj( NewTrajIndxs );
    for i=1:numel(NewTraj)
      trajIndx = NewTrajIndxs(i);
      NewkDists = abs( NewTraj(i) - gridKs );
      NewShortDistIndxs = find( NewkDists < kDistThresh );
      NewShortDists = NewkDists( NewShortDistIndxs );
      NewCVals = interp1( kC, C, NewShortDists, 'linear', 0 );
      NewKVals = fftData( NewShortDistIndxs );
      out(trajIndx) = out(trajIndx) + sum( NewKVals .* NewCVals );
    end
  end

end

