
function out = grid_2D( k, traj, N, varargin )
  % out = grid_2D( k, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )

  preImg = iGridT_2D( k, traj, N, varargin{:} );  % create the image


  % Perform density correction
  nK = numel(k);
  nOut = numel(out);
  G = sparse( nOut, nK );
  for i=1:nOut
    NEED TO FIX THIS
  end
  w = lsqr( G, ones(nK,1) );

  out = preImg .* w;
end

