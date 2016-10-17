
function traj = makeTrajPts( nDim, type, varargin )
  % traj = makeTrajPts( nDim, type, parameters );
  %
  % traj is an N x nDim array, where N is the number of trajectory points
  %
  % Inputs:
  % nDim - the number of dimensions
  % type - the type of trajectory to create
  %        'poissonDisc' - poisson disc sampling
  %        'random' - (default)
  %        'radial' - evenly spaced radial spokes
  %
  % Parameters for trajectories:
  % 'poissonDisc': traj = makeTrajPts( nDim, 'poissonDisc', radius );
  %   radius is the radius of the disc (nominal distance between points)
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
    case 'propeller'
      traj = makeTrajPts_propeller( nDim, varargin{:} );
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


function traj = makeTrajPts_propeller( nDim, nReadout, nLines, dkLine, nAngles )
  if nDim ~= 2, error('Not yet implemented'); end;

  nKperAngle = nReadout * nLines;
  dAngle = pi/nAngles;
  
  dkx = 1 / (nReadout+1);
  kx = ones(nLines,1) * linspace(-0.5,0.5,nReadout);
  kyExtent = dkLine * (nLines-1);
  kyMax = kyExtent/2;
  ky = (-kyMax:dkLine:kyMax)' * ones(1,nReadout);

  traj = zeros(nReadout*nLines*nAngles,2);
  for i=1:nAngles
    thisAngle = (i-1)*dAngle;
    thisKx = cos(thisAngle)*kx(:) - sin(thisAngle)*ky(:);
    thisKy = sin(thisAngle)*kx(:) + cos(thisAngle)*ky(:);

    traj((i-1)*nKperAngle+1:i*nKperAngle,1) = thisKx;
    traj((i-1)*nKperAngle+1:i*nKperAngle,2) = thisKy;
  end
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


