
function traj = makeTrajPts( nDim, type, varargin )
  % traj = makeTrajPts( nDim, type, parameters );
  %
  % traj is an N x nDim array, where N is the number of trajectory points
  %
  % Inputs:
  % nDim - the number of dimensions
  % type - the type of trajectory to create
  %        'random' - (default)
  %        'radial' - evenly spaced radial spokes
  %
  % Parameters for trajectories:
  % 'random': traj = makeTrajPts( nDim, 'random', nTraj );
  %   nTraj is the number of points in the trajectory
  % 'radial': traj = makeTrajPts( nDim, 'radial', nSpokes, nPtsPerSpoke );
  %
  % To view the sampling pattern:
  % plot( kTraj(:,1), kTraj(:,2), 'o', 'MarkerFaceColor', 'k', ...
  %   'MarkerEdgeColor', 'k', 'MarkerSize', 4 );
  % set( gca, 'xTick', [], 'yTick', [] );
  %
  % Written by Nicholas Dwork - Copyright 2016

  if nargin<2, type=''; end;

  switch type
    case 'poissonDisc'
      traj = makeTrajPts_poissonDisc( nDim, varargin{:} );
    case 'spinWarp'
      traj = makeTrajPts_spinWarp( nDim, varargin{:} );
    case 'random'
      traj = makeTrajPts_random(nDim, varargin{:} );
    case 'radial'
      traj = makeTrajPts_radial(nDim, varargin{:} );
    otherwise
      traj = makeTrajPts_random(nDim, varargin{:} );
  end

end


function traj = makeTrajPts_poissonDisc( nDim, r )
  if nDim ~= 2, error('Not yet implemented'); end;
  [kx,ky] = randpd2d(r);
  traj = zeros( numel(kx), 2 );
  traj(:,1) = kx;
  traj(:,2) = ky;
end


function traj = makeTrajPts_spinWarp( nDim, dk )
  k = transpose( -0.5 : dk : 0.5-dk );
  traj = zeros( numel(k)^nDim, nDim );

  if nDim == 1
    traj = k;
  elseif nDim == 2
    [kx,ky] = meshgrid( k, k );
    traj(:,1) = kx(:);
    traj(:,2) = ky(:);
  elseif nDim==3
    [kx,ky,kz] = meshgrid( k, k, k );
    traj(:,1) = kx(:);
    traj(:,2) = ky(:);
    traj(:,3) = kz(:);
  end
end


function traj = makeTrajPts_random( nDim, nTraj )
  traj = rand( nTraj, nDim ) - 0.5;
end


function traj = makeTrajPts_radial( nDim, nSpokes, nPtsPerSpoke )
  if nDim ~= 2, error('Radial trajector just for 2 dimensions'); end;

  thetas = linspace( 0, 2*pi, nSpokes );
  rs = linspace( 0, 0.5, nPtsPerSpoke );

  traj = zeros( nSpokes*nPtsPerSpoke, 2 );
  for i=1:nSpokes
    lowIndx = (i-1)*nPtsPerSpoke+1;
    highIndx = i*nPtsPerSpoke;
    traj(lowIndx:highIndx,1) = rs * cos(thetas(i));
    traj(lowIndx:highIndx,2) = rs * sin(thetas(i));
  end
end


