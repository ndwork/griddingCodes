
function [weights,lsqrFlag,lsqrRes] = makePrecompWeights_2D( ...
  traj, N, varargin )
  % [weights,lsqrFlag,lsqrRes] = makePrecompWeights_2D( traj, N, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Determine the density pre-compensation weights to be used in gridding
  %
  % Inputs:
  %   traj is a Mx2 element array specifying the k-space trajectory.
  %     The first/second column is the ky/kx location.
  %     The units are normalized to [-0.5,0.5)
  %   N is a 2 element array [Ny Nx] representing the number of grid points
  %
  % Optional Inputs:
  %   alpha is the oversampling factor > 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Outputs:
  %   weights - 1D array with density compensation weights
  %
  % Optional Outputs:
  %   lsqrFlag - flag describing results of lsqr optimization (see lsqr
  %     documentation)
  %   lsqrRes - lsqr residual
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  iteration = 0;
  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, nGrid );
      out = iGrid_2D( in, traj, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      disp(['lsqr working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
    else
      out = iGridT_2D( in, traj, nGrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );
      out = out(:);
    end
  end

  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b);
  tolerance = 1d-5;
  maxIter = 1000;
  [weights,lsqrFlag,lsqrRes] = lsqr( @applyA, b(:), tolerance, maxIter );
end

