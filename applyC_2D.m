
function out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx, gridKs, varargin )
  % out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx [, gridKs, type, 'noCirc' ] )
  %
  % Inputs:
  % type:  by default, performs a circular convolution  If type=='noCirc',
  % then it performs a regular(non-circular) convolution.
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultType = [];
  p = inputParser;
  p.addParameter( 'type', defaultType );
  p.parse( varargin{:} );
  type = p.Results.type;
  
  if nargin < 8
    gridKs = size2fftCoordinates( N );
    gridKy=gridKs{1};  gridKx=gridKs{2};
    [gridKx,gridKy] = meshgrid(gridKx,gridKy);
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
  end

  nTraj = size(kTraj,1);
  kws = [ max(kCy), max(kCx) ];
  kDistThreshY = kws(1);
  kDistThreshX = kws(2);
  out = zeros( size(gridKy) );
  for trajIndx=1:nTraj
    distsKy = abs( kTraj(trajIndx,1) - gridKy );
    distsKx = abs( kTraj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThreshY & ...
                           distsKx < kDistThreshX );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    out(shortDistIndxs) = out(shortDistIndxs) + ...
      F(trajIndx) * ( CValsY .* CValsX );
  end

  if ~strcmp( type, 'noCirc' )
    for dim=1:2
      alt = zeros( size(kTraj) );

      for altDir=[-1 1]
        alt(:,dim) = altDir;
        newTraj = kTraj + alt;
        if altDir < 0
          newTrajIndxs = find( newTraj(:,dim) > -0.5-kws(dim) );
        else
          newTrajIndxs = find( newTraj(:,dim) < 0.5+kws(dim) );
        end

        newTraj = newTraj( newTrajIndxs, : );
        for i=1:numel(newTrajIndxs)
          trajIndx = newTrajIndxs(i);
          NewDistsKy = abs( newTraj(i,1) - gridKy );
          NewDistsKx = abs( newTraj(i,2) - gridKx );
          NewShortDistIndxs = find( NewDistsKy < kDistThreshY & ...
                                    NewDistsKx < kDistThreshX );
          NewShortDistsKy = NewDistsKy( NewShortDistIndxs );
          NewShortDistsKx = NewDistsKx( NewShortDistIndxs );
          NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
          NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
          out(NewShortDistIndxs) = out(NewShortDistIndxs) + ...
            F(trajIndx) * ( NewCValsY .* NewCValsX );
        end
      end
    end
  end

end
