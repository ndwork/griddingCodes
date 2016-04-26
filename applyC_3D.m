
function out = applyC_3D( F, traj, N, kws, kCy, kCx, kCz, Cy, Cx, Cz )
  % out = applyC_2D( F, traj, N, kws, kCy, kCx, kCz, Cy, Cx, Cz )
  %
  % Written by Nicholas Dwork - Copyright 2016

  Ny = N(1);  kwy=kws(1);
  Nx = N(2);  kwx=kws(2);
  Nz = N(3);  kwz=kws(3);

  gridKs = size2fftCoordinates([ Ny Nx Nz ]);
  gridKy=gridKs{1};  gridKx=gridKs{2};  gridKz=gridKs{3};

  nTraj = size(traj,1);
  kDistThreshY = 0.5*kwy;
  kDistThreshX = 0.5*kwx;
  kDistThreshZ = 0.5*kwz;
  out = zeros( [Ny Nx Nz] );
  for trajIndx = 1:nTraj
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    distsKz = abs( traj(trajIndx,3) - gridKz );
    shortDistIndxsY = find( distsKy < kDistThreshY );
    shortDistIndxsX = find( distsKx < kDistThreshX );
    shortDistIndxsZ = find( distsKz < kDistThreshZ );
    shortDistsKy = distsKy( shortDistIndxsY );
    shortDistsKx = distsKx( shortDistIndxsX );
    shortDistsKz = distsKz( shortDistIndxsZ );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    CValsZ = interp1( kCz, Cz, shortDistsKz, 'linear', 0 );
    CVals = bsxfun( @times, CValsY*transpose(CValsX), ...
      reshape( CValsZ, [1 1 numel(CValsZ)] ) );
    out(shortDistIndxsY,shortDistIndxsX,shortDistIndxsZ ) = ...
      out(shortDistIndxsY,shortDistIndxsX,shortDistIndxsZ ) + ...
      F(trajIndx) * CVals;
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
        NewShortDistIndxsY = find( NewDistsKy < kDistThreshY );
        NewShortDistIndxsX = find( NewDistsKx < kDistThreshX );
        NewShortDistIndxsZ = find( NewDistsKz < kDistThreshZ );
        NewShortDistsKy = NewDistsKy( NewShortDistIndxsY );
        NewShortDistsKx = NewDistsKx( NewShortDistIndxsX );
        NewShortDistsKz = NewDistsKz( NewShortDistIndxsZ );
        NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
        NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
        NewCValsZ = interp1( kCz, Cz, NewShortDistsKz, 'linear', 0 );
        NewCVals = bsxfun( @times, NewCValsY*transpose(NewCValsX), ...
          reshape( NewCValsZ, [1 1 numel(NewCValsZ)] ) );
        out( NewShortDistIndxsY, NewShortDistIndxsX, NewShortDistIndxsZ ) = ...
          out( NewShortDistIndxsY, NewShortDistIndxsX, NewShortDistIndxsZ ) + ...
          F(trajIndx) * NewCVals;
      end
    end
  end

end
