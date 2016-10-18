
function out = applyC2Grid_2D( F, kTraj, N, kCx, kCy, Cx, Cy )
  % out = applyC2Grid_2D( F, kTraj, N, kCx, kCy, Cx, Cy )
  %
  % Inputs:
  % F - 1D array of Fourier values
  % kTraj - Mx2 array where first/second column represents ky/kx location
  %   of trajectory
  % N - 2 element array specifying grid size [Ny Nx]
  %
  % Written by Nicholas Dwork - Copyright 2016
  % Based on codes written by John Pauly and Ethan Johnson

  kws = [ max(kCy)*N(1), max(kCx)*N(2) ];  % kernel widths
  rkws = round( kws );

  % Convert k-space coordinates to (non-integer) indexes into grid
  kTrajIndxs = zeros( size(kTraj) );
  kTrajIndxs(:,1) = (floor(N(1)/2)+1) + N(1)*kTraj(:,1);
  kTrajIndxs(:,2) = (floor(N(2)/2)+1) + N(2)*kTraj(:,2);
  kCyIndxs = kCy * N(1)/2;
  kCxIndxs = kCx * N(2)/2;

  out = zeros(N);
  for dky = -(rkws(1)-1)/2:(rkws(1)-1)/2
    nearOutKyIndxs = round( kTrajIndxs(:,1) + dky );
    kyDists = abs( nearOutKyIndxs - kTrajIndxs(:,1) );
    cValsKy = interp1( kCyIndxs, Cy, kyDists, 'linear', 0 );
    nearOutKyIndxs = mod( nearOutKyIndxs-1, N(1) ) + 1;

    for dkx = -(rkws(2)-1)/2:(rkws(2)-1)/2
      nearOutKxIndxs = round( kTrajIndxs(:,2) + dkx );
      kxDists = abs( nearOutKxIndxs - kTrajIndxs(:,2) );
      cValsKx = interp1( kCxIndxs, Cx, kxDists, 'linear', 0 );
      nearOutKxIndxs = mod( nearOutKxIndxs-1, N(2) ) + 1;

      out = out + sparse( nearOutKxIndxs, nearOutKyIndxs, ...
        F .* cValsKx .* cValsKy, N(2), N(1) );
    end
  end

end
