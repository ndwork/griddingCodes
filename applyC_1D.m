
function out = applyC_1D( F, traj, N, kw, kC, C )
  % out = applyC_1D( F, traj, N, kw, kC, C )
  %
  % Written by Nicholas Dwork - Copyright 2016

  gridKs = size2fftCoordinates( N );
  nTraj = numel(traj);
  kDistThresh = 0.5*kw;
  out = zeros( N, 1 );
  for trajIndx=1:nTraj
    kDists = abs( traj(trajIndx) - gridKs );
    shortDistIndxs = find( kDists < kDistThresh );
    shortDists = kDists( shortDistIndxs );
    CVals = interp1( kC, C, shortDists, 'linear', 0 );
    out(shortDistIndxs) = out(shortDistIndxs) + F(trajIndx) * CVals;
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
      out(NewShortDistIndxs) = out(NewShortDistIndxs) + ...
        F(trajIndx) * NewCVals;
    end
  end

end
