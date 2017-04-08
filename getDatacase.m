

function [kTraj,fVals,N,psfMask] = getDatacase( datacase, dataDir, varargin )
  % [kTraj,fVals,N,psfMask] = getDatacase( datacase )
  %
  % datacase 0 is a simulation of NickPhantom

  psfMask = [];

  switch datacase

    case 0
      % Simulated data of NickPhantom
      [fVals,phantImg] = makeNickPhantom( kTraj );                               %#ok<ASGLU>
      N = [256, 256];
      psfMask = ones(2*N);

    case 1
      dataFile = 'data_2drad_test.mat';
      load( [dataDir dataFile] );
      ds = 4;  % Downsample factor
      fVals = d(1:end,1:ds:end);
      fVals = fVals(:);

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
      fVals = zeros( prod(sData(1:2)), nCoils ); 
      for cIndx=1:nCoils
        fVals(:,cIndx) = reshape( d(:,:,cIndx), [prod(sData(1:2)) 1] );
      end
      ks = k(:);
      kTraj = [ real(ks) imag(ks) ];

      %N = [ 128 64 ];
      %N = [ 150 74 ];
      %N = [ 140 70 ];
      %N = [ 128 128 ];
      N = [ 150 150 ];
      psfMask = ones(2*N);

    case 3
      dataFile = 'propeller_phantom_2015_11_05.mat';
      load( [dataDir dataFile] );

      ds = 1;  % Downsample factor
      data = data(81:ds:239,1:2:end,1:22,:);
      ktraj = ktraj(81:ds:239,1:2:end,1:22) * 2;
      
      sData = size( data );
      nCoils = sData(4);
      fVals = zeros( prod(sData(1:3)), nCoils );
      for cIndx=1:nCoils
        fVals(:,cIndx) = reshape( data(:,:,:,cIndx), [prod(sData(1:3)) 1] );
      end
      kTraj = [ real(ktraj(:)) imag(ktraj(:)) ];

      N = [ 150 150 ];
      radImg = makeRadialImg( 2*N );
      psfMask = abs(radImg) < min(N);
  end

  %figure; scatter( kTraj(:,1), kTraj(:,2), 8, 'k', 'filled' );
end
