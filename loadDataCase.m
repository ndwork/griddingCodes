
function [kTraj,iGridFVals,N,psfMask] = loadDataCase( datacase, varargin )
  % [kTraj,iGridFVals] = loadDataCase( datacase [, ds ] )
  %
  % ds is the amount to downsample the data

  dataDir = '/Volumes/Seagate2TB/Data/griddingData/';
  %dataDir = '/Users/ndwork/Desktop/griddingData/';

  p = inputParser;
  p.addOptional( 'ds', 1, @isnumeric );
  p.parse(varargin{:});
  ds = p.Results.ds;

  psfMask = [];

  switch datacase
    case 1
      dataFile = 'data_2drad_test.mat';
      load( [dataDir dataFile] );
      iGridFVals = d(1:end,1:ds:end);
      iGridFVals = iGridFVals(:);

      ks = k(1:end,1:ds:end);
      ks = ks(:);
      kTraj = [ real(ks) imag(ks) ];

      N = [ 256 256 ];
      radImg = makeRadialImg( 2*N );
      psfMask = abs(radImg) < min(N);

    case 2
      dataFile = 'data_2dspiralNav.mat';
      load( [dataDir dataFile] );
      sData = size( d );
      nCoils = sData(3);
      iGridFVals = zeros( prod(sData(1:2)), nCoils ); 
      for cIndx=1:nCoils
        iGridFVals(:,cIndx) = reshape( d(:,:,cIndx), [prod(sData(1:2)) 1] );
      end
      ks = k(:);
      kTraj = [ real(ks) imag(ks) ];

      %N = [ 128 64 ];
      N = [ 128 128 ];

    case 3
      dataFile = 'propeller_phantom_2015_11_05.mat';
      load( [dataDir dataFile] );

      data = data(1:3:end,1:5:end,1:4:end,:);
      ktraj = ktraj(1:3:end,1:5:end,1:22);

      sData = size( data );
      nCoils = sData(4);
      iGridFVals = zeros( prod(sData(1:3)), nCoils );
      for cIndx=1:nCoils
        iGridFVals(:,cIndx) = reshape( data(:,:,:,cIndx), [prod(sData(1:3)) 1] );
      end
      kTraj = [ real(ktraj(:)) imag(ktraj(:)) ];

      N = [ 512 512 ];
  end

  %figure; scatter( kTraj(:,1), kTraj(:,2), 8, 'k', 'filled' );
end