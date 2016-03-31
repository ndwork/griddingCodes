

function testModules
  close all; clear; rng(1);

  %% Test iGrid_1D
  nPts = 100;
  rectWidth = 31;
  dataCoords = (0:(nPts-1)) - ceil((nPts-1)/2);
  rect = double( abs(dataCoords) <= rectWidth/2 )';

  %kTraj = rand( nPts, 1 ) - 0.5;
  dk = 1/nPts;
  kTraj = -0.5 : dk : 0.5-dk;
  trueFourierValues = rectWidth .* sinc( rectWidth .* kTraj );
  iGridFourierValues = iGrid_1D( rect, kTraj, 'alpha', 1.5, 'W', 8 );
  error = norm( trueFourierValues(:) - iGridFourierValues(:), 2 ) / ...
          norm( trueFourierValues, 2 );
  disp([ 'iGrid_1D error: ', num2str(error) ]);

  %% Make sure iGrid_1D and iGridT_1D are adjoints
  nXPts = 10;
  nYPts = 5;
  kTraj = rand( nYPts, 1 ) - 0.5;
  %dk = 1/nPts;
  %kTraj = -0.5 : dk : 0.5-dk;
  x = rand( nXPts, 1 );
  y = rand( nYPts, 1 );
  Ax = iGrid_1D( x, kTraj );
  innerProd1 = dotP( Ax, y );
  ATy = iGridT_1D( y, kTraj, nXPts );
  innerProd2 = dotP( x, ATy );
  error = abs( innerProd1 - innerProd2 );
  disp([ 'iGrid/iGridT 1D Adjointness error:  ', num2str(error) ]);

  %% Test iGrid_2D
  % iGrid_2D:  test for values when transforming known object
  % The object is a simple square, which transforms to a 2D sinc
  % The k-space trajector will be randomly selected point
  nKPts = 500;
  nPts = 512;
  squareWidth = 2/100;
  kTraj = rand( nKPts, 2 ) - 0.5;
  imgCoords = (0:(nPts-1)) - ceil(nPts/2);
  [x,y] = meshgrid( imgCoords, imgCoords );
  p = 1/(squareWidth*0.5);
  img = abs(x) < p*0.5 & abs(y) < p*0.5;
  trueKxVals = p .* sinc( p .* kTraj(:,1) );
  trueKyVals = p .* sinc( p .* kTraj(:,2) );
  trueKVals = trueKxVals .* trueKyVals;
  tic;
  kVals = iGrid_2D( img, kTraj );
  iGridTime = toc;
  disp(['iGrid_2D time: ', num2str(iGridTime)]);
  error = norm( trueKVals - kVals, 2 ) / norm( trueKVals, 2 );
  disp(['iGrid_2D error: ', num2str(error)]);

  %% Make sure iGrid_2D and iGridT_2D are adjoints
  sizeX = [ 50, 50 ];
  nY = 20;
  kTraj = rand( nY, 2 ) - 0.5;
  x = rand( sizeX );
  y = rand( nY, 1 );
  Ax = iGrid_2D( x, kTraj );
  innerProd1 = dotP( Ax, y );
  ATy = iGridT_2D( y, kTraj, sizeX );
  innerProd2 = dotP( x, ATy );
  error = abs( innerProd1 - innerProd2 );
  disp([ 'iGrid/iGridT 2D Adjointness error:  ', num2str(error) ]);

end

