
function out = applyCT_2D( F, kTraj, N, kCy, kCx, Cy, Cx, gridKs )
  % out = applyCT_2D( F, kTraj, N, kCy, kCx, Cy, Cx [, gridKs ] )
  %
  % Inputs:
  % fftData - Mx1 array specifying Fourier value at traj(indx,:)
  % traj - Mx2 array where first/second column represents ky/kx location
  %   of trajectory
  % N - 2 element array specifying [Ny Nx] points in grid
  %
  % Written by Nicholas Dwork - Copyright 2016

  if nargin < 8
    gridKs = size2fftCoordinates( N );
    gridKy=gridKs{1};  gridKx=gridKs{2};
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
  end

  kDistThreshY = max(kCy);
  kDistThreshX = max(kCx);

  nTraj = size( kTraj, 1 );
  out = zeros( nTraj, 1 );
  indxs = 1:size(kTraj,1);
  for j=1:numel(gridKy)
    distsKy = abs( gridKy(j) - kTraj(:,1) );
    shortDistsKy = distsKy( distsKy < kDistThreshY );
    if numel(shortDistsKy)==0, continue; end;

    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    subTraj = kTraj( distsKy < kDistThreshY, : );
    subIndxsY = indxs( distsKy < kDistThreshY );

    for i=1:numel(gridKx)
      distsKx = abs( gridKx(i) - subTraj(:,2) );
      shortDistsKx = distsKx( distsKx < kDistThreshX );
      if numel(shortDistsKx)==0, continue; end;

      subCValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
      subCValsY = CValsY( distsKx < kDistThreshX );
      subIndxsYX = subIndxsY( distsKx < kDistThreshX );

      out(subIndxsYX) = out(subIndxsYX) + F(j,i) * subCValsX .* subCValsY;
    end
  end

  kws = [ kDistThreshY, kDistThreshX ];
  onesCol = ones(nTraj,1);
  for dim=1:2
    alt = zeros( size(kTraj) );
    alt(:,dim) = onesCol;

    for altDir=[-1 1]
      NewTraj = kTraj + altDir*alt;
      if altDir < 0
        NewTrajIndxs = find( NewTraj(:,dim) > -0.5-kws(dim) );
      else
        NewTrajIndxs = find( NewTraj(:,dim) < 0.5+kws(dim) );
      end

      NewTraj = NewTraj( NewTrajIndxs, : );
      for j=1:numel(gridKy)
        distsKy = abs( gridKy(j) - NewTraj(:,1) );
        shortDistsKy = distsKy( distsKy < kDistThreshY );
        if numel(shortDistsKy)==0, continue; end;

        CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
        subTraj = NewTraj( distsKy < kDistThreshY, : );
        subIndxsY = NewTrajIndxs( distsKy < kDistThreshY );

        for i=1:numel(gridKx)
          distsKx = abs( gridKx(i) - subTraj(:,2) );
          shortDistsKx = distsKx( distsKx < kDistThreshX );
          if numel(shortDistsKx)==0, continue; end;

          subCValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
          subCValsY = CValsY( distsKx < kDistThreshX );
          subIndxsYX = subIndxsY( distsKx < kDistThreshX );

          out(subIndxsYX) = out(subIndxsYX) + F(j,i) * subCValsX .* subCValsY;
        end
      end

    end
  end

end
