

function testModules
  close all; clear; rng(1);

  %% Test 2D DFT
  nPts = 100;
  rectWidth = 31;
  imgCoords = size2imgCoordinates( [nPts nPts] );
  [x,y] = meshgrid( imgCoords{1}, imgCoords{2} );
  rect2D = double( abs(x) <= rectWidth/2 & abs(y) <= rectWidth/2 )';
  dk = 1/nPts;
  ky = transpose( -0.5 : dk : 0.5-dk );
  kx = transpose( -0.5 : dk : 0.5-dk );
  [kx,ky] = meshgrid( ky, kx );
  kTraj = zeros(nPts*nPts,2);
  kTraj(:,1) = ky(:);
  kTraj(:,2) = kx(:);
  trueFyVals = rectWidth .* sinc( rectWidth .* kTraj(:,1) );
  trueFxVals = rectWidth .* sinc( rectWidth .* kTraj(:,2) );
  trueFVals = trueFyVals .* trueFxVals;
  DFTValues = fftshift( fft2( ifftshift(rect2D) ) );
  error = norm( trueFVals(:) - DFTValues(:), 2 ) / ...
          norm( trueFVals(:), 2 );
  disp([ '2D iDFT error: ', num2str(error) ]);

  %% Test iGrid_2D
  nPts = 100;
  rectWidth = 31;
  imgCoords = size2imgCoordinates( [nPts nPts] );
  [x,y] = meshgrid( imgCoords{1}, imgCoords{2} );
  rect2D = double( abs(x) <= rectWidth/2 & abs(y) <= rectWidth/2 )';
  nTraj = 10000;
  kTraj = rand( nTraj, 2 ) - 0.5;
  trueFyVals = rectWidth .* sinc( rectWidth .* kTraj(:,1) );
  trueFxVals = rectWidth .* sinc( rectWidth .* kTraj(:,2) );
  trueFVals = trueFyVals .* trueFxVals;
  iGridFVals = iGrid_2D( rect2D, kTraj );
  error = norm( trueFVals - iGridFVals, 2 ) / norm( trueFVals, 2 );
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


  %% Test DFT
  nPts = 100;
  rectWidth = 31;
  dataCoords = (0:(nPts-1)) - ceil((nPts-1)/2);
  rect = double( abs(dataCoords) <= rectWidth/2 )';
  dk = 1/nPts;
  kTraj = -0.5 : dk : 0.5-dk;
  trueFourierValues = rectWidth .* sinc( rectWidth .* kTraj );
  DFTValues = fftshift( fft( ifftshift(rect) ) );
  error = norm( trueFourierValues(:) - DFTValues(:), 2 ) / ...
          norm( trueFourierValues, 2 );
  disp([ 'iDFT error: ', num2str(error) ]);

  %% Test grid_1D
  nTraj = 500;
  rectWidth = 31;
  nPts = 100;
  dataCoords = (0:(nPts-1))' - ceil((nPts-1)/2)';
  trueRect = double( abs(dataCoords) <= rectWidth/2 )';
  kTraj = rand( nTraj, 1 ) - 0.5;
  %dk = 1/nTraj;
  %kTraj = (-0.5 : dk : 0.5-dk)';
  F = rectWidth .* sinc( rectWidth .* kTraj );
  weights = makePrecompWeights_1D( kTraj, nPts );
  gridRect = grid_1D( F, kTraj, nPts, weights );
  err = norm( trueRect(:) - gridRect(:), 2 );
  disp([ 'grid_1D error: ', num2str(err) ]);

  %% Test iGrid_1D
  nPts = 100;
  rectWidth = 31;
  dataCoords = (0:(nPts-1)) - ceil((nPts-1)/2);
  rect = double( abs(dataCoords) <= rectWidth/2 )';
  nTraj = 1000;
  kTraj = rand( nTraj, 1 ) - 0.5;
  %dk = 1/nTraj;
  %kTraj = -0.5 : dk : 0.5-dk;
  trueFourierValues = rectWidth .* sinc( rectWidth .* kTraj );
  iGridFourierValues = iGrid_1D( rect, kTraj );
  err = norm( trueFourierValues(:) - iGridFourierValues(:), 2 ) / ...
          norm( trueFourierValues, 2 );
  disp([ 'iGrid_1D error: ', num2str(err) ]);

  % Test grid_1D - roundtrip with iGrid_1D
  weights = makePrecompWeights_1D( kTraj, nPts );
  gridded_1D = grid_1D( iGridFourierValues, kTraj, nPts, weights, ...
    'alpha', 1.5, 'W', 8 );
  err = norm( gridded_1D(:) - rect(:), 2 ) / ...
    norm( rect(:), 2 );
  disp([ 'grid_1D roundtrip error: ', num2str(err) ]);
  %figure; subplot(1,2,1); scatter( kTraj, weights );
  %subplot(1,2,2); plot( rect, 'k' ); hold on; plot( real(gridded_1D),'b' ); plot( imag(gridded_1D),'r' );

  %% Make sure iGrid_1D and iGridT_1D are adjoints
  nXPts = 100;
  nYPts = 50;
  kTraj = rand( nYPts, 1 ) - 0.5;
  %dk = 1/nYPts;
  %kTraj = -0.5 : dk : 0.5-dk;
  x = rand( nXPts, 1 );
  y = rand( nYPts, 1 );
  Ax = iGrid_1D( x, kTraj );
  innerProd1 = dotP( Ax, y );
  ATy = iGridT_1D( y, kTraj, nXPts );
  innerProd2 = dotP( x, ATy );
  error = abs( innerProd1 - innerProd2 );
  disp([ 'iGrid/iGridT 1D Adjointness error:  ', num2str(error) ]);



%   %% Test grid_2D
%   %  Test with Fourier Transform of known 2D values
  imgFile = '/Applications/MATLAB_R2013a_Student.app/toolbox/images/imdemos/moon.tif';
%   img = double( imread( imgFile ) );
%   img = img(1:5:end,1:5:end);
%   fftImg = fftshift( fft2( ifftshift( img ) ) );
%   ks = size2fftCoordinates( size( fftImg ) );
%   ky=ks{1};  kx=ks{2};
%   [kx,ky] = meshgrid( kx, ky );
%   traj = zeros( numel(img), 2 );
%   traj(:,1) = ky(:);  traj(:,2) = kx(:);
%   out = grid_2D( fftImg(:), traj, size(img), 'alpha', 1.5, 'W', 10 );
%   figure; imshow(imresize(img,3,'nearest'),[0 255]); title('origImg');
%   figure; imshow(imresize(out,3,'nearest'),[0 255]); title('out');
%   error = norm( out(:) - img(:), 2 ) / norm( img(:), 2 );
%   disp([ 'iGridT_2D error: ', num2str(error) ]);






  %% Test iGrid_3D
  % iGrid_3D:  test for values when transforming known object
  % The object is a simple square, which transforms to a 3D sinc
  % The k-space trajector will be randomly selected point
  nKPts = 500;
  nPts = 512;
  squareWidth = 2/100;
  kTraj = rand( nKPts, 3 ) - 0.5;
  imgCoords = (0:(nPts-1)) - ceil(nPts/2);
  [x,y,z] = meshgrid( imgCoords, imgCoords, imgCoords );
  p = 1/(squareWidth*0.5);
  img = abs(x) < p*0.5 & abs(y) < p*0.5 & abs(z) < p*0.5;
  trueKxVals = p .* sinc( p .* kTraj(:,1) );
  trueKyVals = p .* sinc( p .* kTraj(:,2) );
  trueKzVals = p .* sinc( p .* kTraj(:,3) );
  trueKVals = trueKxVals .* trueKyVals .* trueKzVals;
  tic;
  kVals = iGrid_3D( img, kTraj );
  iGrid_3D_time = toc;
  error = norm( trueKVals - kVals, 2 ) / norm( trueKVals, 2 );
  disp(['iGrid_3D time taken: ', num2str(iGrid_3D_time)]);
  disp(['iGrid_3D error: ', num2str(error)]);

  %% Make sure iGrid_3D and iGridT_3D are adjoints
  sizeX = [ 50, 50, 50 ];
  nY = 20;
  kTraj = rand( nY, 3 ) - 0.5;
  x = rand( sizeX );
  y = rand( nY, 1 );
  Ax = iGrid_3D( x, kTraj );
  innerProd1 = dotP( Ax, y );
  ATy = iGridT_3D( y, kTraj, sizeX );
  innerProd2 = dotP( x, ATy );
  error = abs( innerProd1 - innerProd2 );
  disp([ 'iGrid/iGridT 3D Adjointness error:  ', num2str(error) ]);
end

