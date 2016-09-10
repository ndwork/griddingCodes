
function out = applyC_3D( F, traj, N, kCy, kCx, kCz, Cy, Cx, Cz, gridKs )
  % out = applyC_3D( F, traj, N, kCy, kCx, kCz, Cy, Cx, Cz, gridKs )
  %
  % Written by Nicholas Dwork - Copyright 2016

  if nargin < 10
    gridKs = size2fftCoordinates( N );
    gridKy=gridKs{1};  gridKx=gridKs{2};  gridKz=gridKs{3};
    [gridKx,gridKy,gridKz] = meshgrid(gridKx,gridKy,gridKz);
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
    gridKz = gridKs(:,3);
  end

  nTraj = size( traj, 1 );
  kws = [ max(kCy), max(kCx), max(kCz) ];
  kDistThreshY = kws(1);
  kDistThreshX = kws(2);
  kDistThreshZ = kws(3);
  tmp = cell(nTraj,1);
  parfor trajIndx=1:nTraj
    disp(['applyC_3D: Working on traj ', num2str(trajIndx), ...
      ' of ', num2str(nTraj)]);
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    distsKz = abs( traj(trajIndx,3) - gridKz );
    shortDistIndxs = find( distsKy < kDistThreshY & ...
                           distsKz < kDistThreshX & ...
                           distsKx < kDistThreshZ );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    shortDistsKz = distsKz( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    CValsZ = interp1( kCz, Cz, shortDistsKz, 'linear', 0 );

    tmp{trajIndx} = struct( ...
      'indxs', shortDistIndxs, ...
      'values', F(trajIndx) * ( CValsY .* CValsX .* CValsZ ) ...
    );
  end

  out = zeros( size(gridKy) );
  for trajIndx=1:nTraj
    out( tmp{trajIndx}.indxs ) = out( tmp{trajIndx}.indxs ) + ...
      tmp{trajIndx}.values;
  end

  onesCol = ones(nTraj,1);
  for dim=1:3
    alt = zeros( size(traj) );
    alt(:,dim) = onesCol;

    for altDir=[-1 1]
      NewTraj = traj + altDir*alt;
      if altDir < 0
        NewTrajIndxs = find( NewTraj(:,dim) > -0.5-kws(dim) );
      else
        NewTrajIndxs = find( NewTraj(:,dim) < 0.5+kws(dim) );
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
        out(NewShortDistIndxs) = out(NewShortDistIndxs) + ...
          F(trajIndx) * ( NewCValsY .* NewCValsX .* NewCValsZ );
      end
    end
  end

end
