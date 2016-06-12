
function areas = voronoidens(kTraj)
%
% function areas = voronoidens(kTraj);
%
% output: areas of cells for each point 
%           (if point doesn't have neighbors the area is NaN)

% uncomment these to plot voronoi diagram
%[vx, vy] = voronoi(kx,ky);
%plot(kx,ky,'r.',vx,vy,'b-'); axis equal

[V,C] = voronoin(kTraj);  % returns vertices and cells of voronoi diagram
nK = size( kTraj, 1 );
areas = zeros( nK, 1 );
for j = 1:nK
  x = V(C{j},1); y = V(C{j},2); lxy = length(x);
  A = abs( sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:)) ) );
  areas(j) = A;
end

areas( ~isfinite(areas) ) = 0;
