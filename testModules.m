

function testModules
  close all; clear; rng(3);

%   %% Test DFT
%   nPts = 100;
%   rectWidth = 31;
%   dataCoords = (0:(nPts-1)) - ceil((nPts-1)/2);
%   rect = double( abs(dataCoords) <= rectWidth/2 )';
%   dk = 1/nPts;
%   kTraj = -0.5 : dk : 0.5-dk;
%   trueFourierValues = rectWidth .* sinc( rectWidth .* kTraj );
%   DFTValues = fftshift( fft( ifftshift(rect) ) );
%   error = norm( trueFourierValues(:) - DFTValues(:), 2 ) / ...
%           norm( trueFourierValues, 2 );
%   disp([ 'iDFT error: ', num2str(error) ]);
% 
% 
%   %% Test iGrid_1D
%   nPts = 100;
%   nTraj = 500;
%   rectWidth = 31;
%   dataCoords = (0:(nPts-1)) - ceil((nPts-1)/2);
%   rect = double( abs(dataCoords) <= rectWidth/2 )';
%   kTraj = rand( nTraj, 1 ) - 0.5;
%   %dk = 1/nTraj;
%   %kTraj = -0.5 : dk : 0.5-dk;
%   trueFourierValues = rectWidth .* sinc( rectWidth .* kTraj );
%   iGridFourierValues = iGrid_1D( rect, kTraj );
%   err = norm( trueFourierValues(:) - iGridFourierValues(:), 2 ) / ...
%           norm( trueFourierValues, 2 );
%   disp([ 'iGrid_1D error: ', num2str(err) ]);
% 
%   % Test grid_1D - roundtrip with iGrid_1D
%   weights = makePrecompWeights_1D( kTraj, nPts );
%   gridded_1D = grid_1D( iGridFourierValues, kTraj, nPts, weights );
%   err = norm( gridded_1D(:) - rect(:), 2 ) / ...
%     norm( rect(:), 2 );
%   disp([ 'grid_1D roundtrip error: ', num2str(err) ]);
%   %figure; subplot(1,2,1); scatter( kTraj, weights );
%   %subplot(1,2,2); plot( rect, 'k' ); hold on; plot( real(gridded_1D),'b' ); plot( imag(gridded_1D),'r' );
% 
%   %% Test grid_1D
%   nTraj = 500;
%   rectWidth = 31;
%   nPts = 100;
%   dataCoords = (0:(nPts-1))' - ceil((nPts-1)/2)';
%   trueRect = double( abs(dataCoords) <= rectWidth/2 )';
%   kTraj = rand( nTraj, 1 ) - 0.5;
%   %dk = 1/nTraj;
%   %kTraj = (-0.5 : dk : 0.5-dk)';
%   F = rectWidth .* sinc( rectWidth .* kTraj );
%   weights = makePrecompWeights_1D( kTraj, nPts );
%   gridRect = grid_1D( F, kTraj, nPts, weights );
%   err = norm( trueRect(:) - gridRect(:), 2 );
%   disp([ 'grid_1D error: ', num2str(err) ]);
% 
% 
%   %% Make sure grid_1D and gridT_1D are adjoints
%   nXPts = 50;
%   nYPts = 100;
%   kTraj = rand( nXPts, 1 ) - 0.5;
%   weights = rand( nXPts, 1 );
%   %dk = 1/nYPts;
%   %kTraj = -0.5 : dk : 0.5-dk;
%   x = rand( nXPts, 1 );
%   y = rand( nYPts, 1 );
%   Ax = grid_1D( x, kTraj, nYPts, weights );
%   ATy = gridT_1D( y, kTraj, nYPts, weights );
%   innerProd1 = dotP( Ax, y );
%   innerProd2 = dotP( x, ATy );
%   error = abs( innerProd1 - innerProd2 );
%   disp([ 'grid/gridT 1D Adjointness error:  ', num2str(error) ]);  
% 
% 
%   %% Make sure iGrid_1D and iGridT_1D are adjoints
%   nXPts = 100;
%   nYPts = 50;
%   kTraj = rand( nYPts, 1 ) - 0.5;
%   %dk = 1/nYPts;
%   %kTraj = -0.5 : dk : 0.5-dk;
%   x = rand( nXPts, 1 );
%   y = rand( nYPts, 1 );
%   Ax = iGrid_1D( x, kTraj );
%   innerProd1 = dotP( Ax, y );
%   ATy = iGridT_1D( y, kTraj, nXPts );
%   innerProd2 = dotP( x, ATy );
%   error = abs( innerProd1 - innerProd2 );
%   disp([ 'iGrid/iGridT 1D Adjointness error:  ', num2str(error) ]);
% 
% 
%   %% Test 2D DFT
%   nPts = 100;
%   rectWidth = 31;
%   imgCoords = size2imgCoordinates( [nPts nPts] );
%   [x,y] = meshgrid( imgCoords{1}, imgCoords{2} );
%   rect2D = double( abs(x) <= rectWidth/2 & abs(y) <= rectWidth/2 )';
%   dk = 1/nPts;
%   ky = transpose( -0.5 : dk : 0.5-dk );
%   kx = transpose( -0.5 : dk : 0.5-dk );
%   [kx,ky] = meshgrid( ky, kx );
%   kTraj = zeros(nPts*nPts,2);
%   kTraj(:,1) = ky(:);
%   kTraj(:,2) = kx(:);
%   trueFyVals = rectWidth .* sinc( rectWidth .* kTraj(:,1) );
%   trueFxVals = rectWidth .* sinc( rectWidth .* kTraj(:,2) );
%   trueFVals = trueFyVals .* trueFxVals;
%   DFTValues = fftshift( fft2( ifftshift(rect2D) ) );
%   error = norm( trueFVals(:) - DFTValues(:), 2 ) / ...
%           norm( trueFVals(:), 2 );
%   disp([ '2D iDFT error: ', num2str(error) ]);
% 
% 
%   %% Make sure grid_1D and gridT_1D are adjoints
%   nXPts = 50;
%   nYPts = 100;
%   kTraj = rand( nXPts, 1 ) - 0.5;
%   weights = rand( nXPts, 1 );
%   %dk = 1/nYPts;
%   %kTraj = -0.5 : dk : 0.5-dk;
%   x = rand( nXPts, 1 );
%   y = rand( nYPts, 1 );
%   Ax = grid_1D( x, kTraj, nYPts, weights );
%   ATy = gridT_1D( y, kTraj, nYPts, weights );
%   innerProd1 = dotP( Ax, y );
%   innerProd2 = dotP( x, ATy );
%   error = abs( innerProd1 - innerProd2 );
%   disp([ 'grid/gridT 1D Adjointness error:  ', num2str(error) ]);
% 
% 
%   %% Make sure iGrid_2D and iGridT_2D are adjoints
%   sizeX = [ 50, 50 ];
%   nY = 20;
%   kTraj = rand( nY, 2 ) - 0.5;
%   x = rand( sizeX );
%   y = rand( nY, 1 );
%   Ax = iGrid_2D( x, kTraj );
%   innerProd1 = dotP( Ax, y );
%   ATy = iGridT_2D( y, kTraj, sizeX );
%   innerProd2 = dotP( x, ATy );
%   error = abs( innerProd1 - innerProd2 );
%   disp([ 'iGrid/iGridT 2D Adjointness error:  ', num2str(error) ]);
% 
% 
%   %% Make sure grid_2D and gridT_2D are adjoints
%   nX = 100;
%   sizeY = [ 50, 50 ];
%   kTraj = rand( nX, 2 ) - 0.5;
%   weights = rand( nX, 1 );
%   x = rand( nX, 1 );
%   y = rand( sizeY );
%   Ax = grid_2D( x, kTraj, sizeY, weights );
%   ATy = gridT_2D( y, kTraj, sizeY, weights );
%   innerProd1 = dotP( Ax, y );
%   innerProd2 = dotP( x, ATy );
%   error = abs( innerProd1 - innerProd2 );
%   disp([ 'iGrid/iGridT 2D Adjointness error:  ', num2str(error) ]);


  %% Test iGrid_2D
  nPts = 64;
  rectWidth = nPts/4 + 1;
  N = [nPts nPts];
  imgCoords = size2imgCoordinates( N );
  [x,y] = meshgrid( imgCoords{1}, imgCoords{2} );
  rect2D = double( abs(x) <= rectWidth/2 & abs(y) <= rectWidth/2 )';
  N = size( rect2D );

  %kTraj = makeTrajPts( 2, 'random', nPts*nPts );
  kTraj = makeTrajPts( 2, 'radial', 360, 50 );
  %kTraj = makeTrajPts( 2, 'spinWarp', 1/nPts );
  %kTraj = makeTrajPts( 2, 'poissonDisc', 1/nPts );
  %plot( kTraj(:,1), kTraj(:,2), 'o', 'MarkerFaceColor', 'k', ...
  %  'MarkerEdgeColor', 'k', 'MarkerSize', 5 ); set( gca, 'xTick', [] ); ...
  %  set( gca, 'yTick', [] );

  ros = makeTrajROS( 'radial', 2*N );
  
  trueFyVals = rectWidth .* sinc( rectWidth .* kTraj(:,1) );
  trueFxVals = rectWidth .* sinc( rectWidth .* kTraj(:,2) );
  trueFVals = trueFyVals .* trueFxVals;
  tic;
  iGridFVals = iGrid_2D( rect2D, kTraj );
  iGrid_2D_time = toc;  %#ok<NASGU>
  error = norm( trueFVals - iGridFVals, 2 ) / norm( trueFVals, 2 );
  disp(['iGrid_2D error: ', num2str(error)]);

  % Test roundtrip error with gridding
  %weights_cls = makePrecompWeights_2D_v2( kTraj, N, ros, 'alg', 'cls' );
  %weights_fp = makePrecompWeights_2D_v2( kTraj, N, 'alg', 'fp' );
  %weights_regCLS = makePrecompWeights_2D_v2( kTraj, N, 'alg', 'regCLS' );
  %weights_regWLS = makePrecompWeights_2D_v2( kTraj, N, 'alg', 'regWLS' );
  weights_wls = makePrecompWeights_2D_v2( kTraj, N, 'alg', 'wls' );
  %weights = makePrecompWeights_2D( kTraj, N, 'alg', 'ls' );
  %weights = makePrecompWeights_2D( kTraj, N, 'alg', 'rls' );
  %weights = makePrecompWeights_2D( kTraj, N, 'alg', 'vor' );
%load( 'weights_64x64_rand.mat' );
  tic;
  gridded_2D = grid_2D( iGridFVals, kTraj, N, weights );
  grid_2D_time = toc;  %#ok<NASGU>
  err = norm( gridded_2D(:) - rect2D(:), 2 ) / ...
    norm( rect2D(:), 2 );
  disp([ 'grid_2D roundtrip error: ', num2str(err) ]);


%   %% Make sure iGrid_2D and iGridT_2D are adjoints
%   sizeX = [ 50, 50 ];
%   nY = 20;
%   kTraj = rand( nY, 2 ) - 0.5;
%   x = rand( sizeX );
%   y = rand( nY, 1 );
%   Ax = iGrid_2D( x, kTraj );
%   innerProd1 = dotP( Ax, y );
%   ATy = iGridT_2D( y, kTraj, sizeX );
%   innerProd2 = dotP( x, ATy );
%   error = abs( innerProd1 - innerProd2 );
%   disp([ 'iGrid/iGridT 2D Adjointness error:  ', num2str(error) ]);
% 
% 
%   %% Test iGrid_3D
%   % iGrid_3D:  test for values when transforming known object
%   % The object is a simple square, which transforms to a 3D sinc
%   % The k-space trajector will be randomly selected point
%   nKPts = 500;
%   nPts = 64;
%   squareWidth = 2/100;
%   kTraj = rand( nKPts, 3 ) - 0.5;
%   rectWidth = 1/(squareWidth*0.5);
%   trueKxVals = rectWidth .* sinc( rectWidth .* kTraj(:,1) );
%   trueKyVals = rectWidth .* sinc( rectWidth .* kTraj(:,2) );
%   trueKzVals = rectWidth .* sinc( rectWidth .* kTraj(:,3) );
%   trueFVals = trueKxVals .* trueKyVals .* trueKzVals;
%   imgCoords = (0:(nPts-1)) - ceil(nPts/2);
%   [x,y,z] = meshgrid( imgCoords, imgCoords, imgCoords );
%   rect3D = double( abs(x) < rectWidth*0.5 & abs(y) < rectWidth*0.5 & ...
%     abs(z) < rectWidth*0.5 );
%   N = size( rect3D );
%   tic;
%   iGridFVals = iGrid_3D( rect3D, kTraj );
%   iGrid_3D_time = toc;
%   error = norm( trueFVals - iGridFVals, 2 ) / norm( trueFVals, 2 );
%   disp(['iGrid_3D time taken: ', num2str(iGrid_3D_time)]);
%   disp(['iGrid_3D error: ', num2str(error)]);
% 
%   % Test roundtrip error with gridding
%   %weights = makePrecompWeights_3D( kTraj, N );
% load 'junk3D.mat';
%   tic;
%   gridded_3D = grid_3D( iGridFVals, kTraj, N, weights );
%   grid_3D_time = toc;
%   err = norm( gridded_3D(:) - rect3D(:), 2 ) / ...
%     norm( rect3D(:), 2 );
%   disp(['grid_3D time taken: ', num2str(grid_3D_time)]);
%   disp([ 'grid_3D roundtrip error: ', num2str(err) ]);
% 
% 
%   %% Make sure iGrid_3D and iGridT_3D are adjoints
%   sizeX = [ 50, 50, 50 ];
%   nY = 20;
%   kTraj = rand( nY, 3 ) - 0.5;
%   x = rand( sizeX );
%   y = rand( nY, 1 );
%   Ax = iGrid_3D( x, kTraj );
%   innerProd1 = dotP( Ax, y );
%   ATy = iGridT_3D( y, kTraj, sizeX );
%   innerProd2 = dotP( x, ATy );
%   error = abs( innerProd1 - innerProd2 );
%   disp([ 'iGrid/iGridT 3D Adjointness error:  ', num2str(error) ]);


end

