

function testModules
  close all; clear; rng(1);


  %% Test iGridT_1D
  imgFile = '/Applications/MATLAB_R2013a_Student.app/toolbox/images/imdemos/moon.tif';
  img = double( imread( imgFile ) );
  img = img(1:5:end,1:5:end);
  data = img(50,:);
  data = transpose( data ./ max(data) );
  nData = numel(data);
  fftData = fftshift( fft( ifftshift( data ) ) );
  kTraj = size2fftCoordinates( numel(fftData) );
  out = iGridT( fftData, kTraj, nData, 'alpha', 1.5, 'W', 8 );
  error = norm( out - data, 2 ) / norm( data, 2 );
  disp([ 'iGridT_1D error: ', num2str(error) ]);


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
  %kTraj = rand( nPts, 1 ) - 0.5;
  dk = 1/nPts;
  kTraj = -0.5 : dk : 0.5-dk;
  x = rand( nPts, 1 );
  y = rand( nPts, 1 );
  Ax = iGrid_1D( x, kTraj );
  innerProd1 = sum( Ax .* y );
  ATy = iGridT_1D( y, kTraj, nPts );
  innerProd2 = sum( x .* ATy );
  error = abs( innerProd1 - innerProd2 );
  disp([ 'iGrid/iGridT 1D Adjointness error:  ', num2str(error) ]);



%   %% Test iGridT_2D
%   %  Test with Fourier Transform of known 2D values
%   imgFile = '/Applications/MATLAB_R2013a_Student.app/toolbox/images/imdemos/moon.tif';
%   img = double( imread( imgFile ) );
%   img = img(1:5:end,1:5:end);
%   fftImg = fftshift( fft2( ifftshift( img ) ) );
%   ks = size2fftCoordinates( size( fftImg ) );
%   ky=ks{1};  kx=ks{2};
%   [kx,ky] = meshgrid( kx, ky );
%   traj = zeros( numel(img), 2 );
%   traj(:,1) = ky(:);
%   traj(:,2) = kx(:);
%   out = iGridT( fftImg(:), traj, size(img), 'alpha', 1.5, 'W', 8 );
%   figure; imshow(imresize(img,3,'nearest'),[0 25]); title('origImg');
%   figure; imshow(imresize(out,3,'nearest'),[0 25]); title('out');
%   error = norm( out(:) - img(:), 2 ) / norm( img(:), 2 );
%   disp([ 'iGridT_2D error: ', num2str(error) ]);






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
  disp(['Time taken for inverse gridding: ', num2str(iGridTime)]);
  error = norm( trueKVals - kVals, 2 ) / norm( trueKVals, 2 );

  tic;
  recon = iGridT_2D( kVals, kTraj, size(img) );
  iGridTTime = toc;
  disp(['Time taken for inverse gridding transpose: ', num2str(iGridTTime)]);


  innerProd1 = sum( conj(kVals) .* kVals );
  innerProd2 = sum( img(:) .* recon(:) );
  
return




  nTheta = 12;
  dTheta = 2*pi / nTheta;
  nPtsPerTheta = 10;
  nTraj = 1 + nPtsPerTheta + nTheta;
  traj = zeros( nTraj, 2 );
  dk = 0.5 / (nTraj+1);
  ks = dk:dk:0.5-dk;
  nKs = numel(ks);
  for i = 0:nTheta
    trajIndx = i*nKs+1;
    kxs = ks * cos( i*dTheta );
    kys = ks * sin( i*dTheta );
    traj(trajIndx+1:trajIndx+nKs,1) = kxs;
    traj(trajIndx+1:trajIndx+nKs,2) = kys;
  end

  % Make sure that applyET_2D is the adjoint of applyE_2D
  N = 512;
  x = phantom(512);
  y = rand(300,1);
  tic
  Ex = iGrid_2D( x, traj );
  ETime= toc;
  disp(['iGrid Time: ', num2str(ETime)]);
  innerProd1 = Ex' * y;

  tic
  ETy = iGridT_2D( y, traj, N );
  ETTime = toc;
  disp(['iGridT Time: ', num2str(ETTime)]);
  innerProd2 = sum( x(:) .* ETy(:) );

  innerProdDiff = abs( innerProd1 - innerProd2 );
  disp(['Inner Product Difference: ', num2str( innerProdDiff ) ]);
  

end

