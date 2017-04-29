
function areas = voronoidens(kTraj)
  %
  % function areas = voronoidens(kTraj);
  %
  % output: areas of cells for each point 
  %           (if point doesn't have neighbors the area is NaN)
  %
  % uncomment these to plot voronoi diagram
  %[vx, vy] = voronoi(kx,ky);
  %plot(kx,ky,'r.',vx,vy,'b-'); axis equal

  [uniqTraj,~,u2k] = unique( kTraj, 'rows' );

  [V,C] = voronoin(uniqTraj);  % returns vertices and cells of voronoi diagram
  nK = size( kTraj, 1 );
  areas = zeros( nK, 1 );
  for j = 1:size( uniqTraj, 1 )
    if numel(C{j}) == 0, continue; end

    x = V(C{j},1); y = V(C{j},2); lxy = length(x);
    A = abs( sum( ...
      0.5 * ( x([2:lxy 1]) - x(:) ) .* ( y([2:lxy 1]) + y(:)) ) ...
    );

    areas( u2k == j ) = A / sum( u2k == j );
  end

  areas( ~isfinite(areas) ) = max( areas( isfinite(areas) ) );

  plotThis = 0;
  if plotThis > 0
    kx = kTraj(:,2);
    ky = kTraj(:,1);
    [vx, vy] = voronoi(kx,ky);
    plot(kx,ky,'r.',vx,vy,'b-');
    axis([-0.5 0.5 -0.5 0.5], 'equal');
  end

end