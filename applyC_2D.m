
function out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx, gridKs )
  % out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx [, gridKs ] )
  %
  % Written by Nicholas Dwork - Copyright 2016

  if nargin < 8
    gridKs = size2fftCoordinates( N );
    gridKy=gridKs{1};  gridKx=gridKs{2};
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
  end

  nTraj = size(kTraj,1);
  kDistThreshY = max(kCy);
  kDistThreshX = max(kCx);

  tmp10 = zeros(nTraj,2);  tmp10(:,1)=1;
  tmp01 = zeros(nTraj,2);  tmp01(:,2)=1;
  traj10_p = kTraj + tmp10;
  F10_p = F( traj10_p(:,1) < 0.5 + kDistThreshY );
  traj10_p = traj10_p( traj10_p(:,1) < 0.5 + kDistThreshY, : );
  traj10_n = kTraj - tmp10;
  F10_n = F( traj10_n(:,1) > -0.5 - kDistThreshY );
  traj10_n = traj10_n( traj10_n(:,1) > -0.5 - kDistThreshY, : );
  traj01_p = kTraj + tmp01;
  F01_p = F( traj01_p(:,2) < 0.5 + kDistThreshX );
  traj01_p = traj01_p( traj01_p(:,2) < 0.5 + kDistThreshX, : );
  traj01_n = kTraj - tmp01;
  F01_n = F( traj01_n(:,2) > -0.5 - kDistThreshX );
  traj01_n = traj01_n( traj01_n(:,2) > -0.5 - kDistThreshX, : );

  allTraj = [ kTraj; traj10_p; traj10_n; traj01_p; traj01_n; ];
  allF = [ F; F10_p; F10_n; F01_p; F01_n; ];

  out = zeros( numel(gridKy), numel(gridKx) );
  for j=1:numel(gridKy)
    distsKy = abs( gridKy(j) - allTraj(:,1) );
    subDists = distsKy( distsKy < kDistThreshY );
    subFy = allF( distsKy < kDistThreshY );
    subTraj = allTraj( distsKy < kDistThreshY, : );
    CValsY = interp1( kCy, Cy, subDists, 'linear', 0 );

    for i=1:numel(gridKx)
      distsKx = abs( gridKx(i) - subTraj(:,2) );
      shortDistsKx = distsKx( distsKx < kDistThreshX );
      CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );

      out(j,i) = sum( subFy(distsKx < kDistThreshX) .* ...
        CValsY(distsKx < kDistThreshX) .* CValsX );
    end
  end

end
