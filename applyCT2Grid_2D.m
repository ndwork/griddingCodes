
function out = applyCT2Grid_2D( fftData, kTraj, N, kCy, kCx, Cy, Cx )
  % out = applyCT2Grid_2D( fftData, kTraj, N, kCy, kCx, Cy, Cx )
  %
  % Inputs:
  % fftData - Mx1 array specifying Fourier value at traj(indx,:)
  % kTraj - Mx2 array where first/second column represents ky/kx location
  %   of trajectory
  % N - 2 element array specifying [Ny Nx] points in grid
  %
  % Written by Nicholas Dwork - Copyright 2016

  kws = [ max(kCy)*N(1), max(kCx)*N(2) ];  % kernel widths
  rkws = round( kws );
  
  % Convert k-space coordinates to (non-integer) indexes into grid
  kTrajIndxs = zeros( size(kTraj) );
  kTrajIndxs(:,1) = (floor(N(1)/2)+1) + N(1)*kTraj(:,1);
  kTrajIndxs(:,2) = (floor(N(2)/2)+1) + N(2)*kTraj(:,2);
  kCyIndxs = kCy * N(1)/2;
  kCxIndxs = kCx * N(2)/2;

  out = zeros( size(kTraj,1), 1 );
  for dkx = -rkws(2):rkws(2)
    nearOutKxIndxs = round( kTrajIndxs(:,1) + dkx );
    kxDists = abs( nearOutKxIndxs - kTrajIndxs(:,1) );
    cValsKx = interp1( kCxIndxs, Cx, kxDists, 'linear', 0 );
    nearOutKxIndxs = mod( nearOutKxIndxs-1, N(1) ) + 1;

    for dky = -(rkws(1)-1):(rkws(1)-1)
      nearOutKyIndxs = round( kTrajIndxs(:,2) + dky );
      kyDists = abs( nearOutKyIndxs - kTrajIndxs(:,2) );
      cValsKy = interp1( kCyIndxs, Cy, kyDists, 'linear', 0 );
      nearOutKyIndxs = mod( nearOutKyIndxs-1, N(2) ) + 1;

      indxs1D = nearOutKyIndxs + (nearOutKxIndxs-1)*N(1);
      out = out + fftData(indxs1D) .* cValsKx .* cValsKy;
    end
  end

end
