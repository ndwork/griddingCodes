
function out = applyCT_3D( fftData, traj, N, ...
  kCy, kCx, kCz, Cy, Cx, Cz [, gridKs] )
  % out = applyCT_3D( fftData, traj, N, kCy, kCx, kCz, Cy, Cx, Cz )
  %
  % Inputs:
  % fftData - Mx1 array specifying Fourier value at traj(indx,:)
  % traj - Mx3 array where first/second/third column represents ky/kx/kz
  %   location of trajectory
  % N - 3 element array specifying [Ny Nx Nz] points in grid
  %
  % Written by Nicholas Dwork - Copyright 2016

  if nargin < 10
    gridKs = size2fftCoordinates( N );
    gridKy=gridKs{1};  gridKx=gridKs{2};  gridKz=gridKs{3};
    [gridKx,gridKy,gridKz] = meshgrid(gridKx,gridKy,gridKz);
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
    gridKz = gridKs(:,2);
  end

  kws = [ max(kCy), max(kCx), max(kCz) ];
  kDistThreshY = kws(1);
  kDistThreshX = kws(2);
  kDistThreshZ = kws(3);

  nTraj = size( traj, 1 );
  out = zeros( nTraj, 1 );
  parfor trajIndx = 1:nTraj
    kDistsY = abs( traj(trajIndx,1) - gridKy );
    kDistsX = abs( traj(trajIndx,2) - gridKx );
    kDistsZ = abs( traj(trajIndx,3) - gridKz );
    shortDistIndxs = find( kDistsY < kDistThreshY & ...
                           kDistsX < kDistThreshX & ...
                           kDistsZ < kDistThreshZ );
    CValsY = interp1( kCy, Cy, shortDistsIndxs, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsIndxs, 'linear', 0 );
    CValsZ = interp1( kCz, Cz, shortDistsIndxs, 'linear', 0 );
    kVals = fftData( shortDistIndxs );
    out( trajIndx ) = sum( kVals .* CValsY .* CValsX .* CValsZ );
  end

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
        out(trajIndx) = out(trajIndx) + sum( NewKVals .* ...
          NewCValsY .* NewCValsX .* NewCValsZ );
      end
    end
  end

end
