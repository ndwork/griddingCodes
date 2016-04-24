
function out = applyCT_2D( fftData, traj, N, kws, kCy, kCx, Cy, Cx )
  % out = applyCT_2D( fftData, traj, N, kws, kCy, kCx, Cy, Cx )
  %
  % Inputs:
  % fftData - Mx1 array specifying Fourier value at traj(indx,:)
  % traj - Mx2 array where first/second column represents ky/kx location
  %   of trajectory
  % N - 2 element array specifying [Ny Nx] points in grid
  %
  % Written by Nicholas Dwork - Copyright 2016

  gridKs = size2fftCoordinates( N );
  gridKy=gridKs{1};  gridKx=gridKs{2};
  [gridKx,gridKy] = meshgrid( gridKx, gridKy );

  kDistThreshY = 0.5*kws(1);
  kDistThreshX = 0.5*kws(2);

  nTraj = size( traj, 1 );
  out = zeros( nTraj, 1 );
  for trajIndx = 1:nTraj
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThreshY & distsKx < kDistThreshX );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    out( trajIndx ) = sum( kVals .* CValsY .* CValsX );
  end

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
        out(trajIndx) = out(trajIndx) + sum( NewKVals .* NewCValsY .* NewCValsX );
      end
    end
  end

end
