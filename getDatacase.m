

function [kTraj,fVals,N,psfMask,phantImg] = getDatacase( datacase, dataDir )
  % [kTraj,fVals,N,psfMask] = getDatacase( datacase )
  %
  % datacase 0 is a simulation of NickPhantom with radial trajectory
  %
  % Written by Nicholas Dwork - Copyright 2017

  phantImg = [];

  switch datacase

    case 0
      % Simulated data of NickPhantom with propeller trajectory
      nReadout = 200;
      nLines = 11;
      dkLine = 0.1;
      nAngles = 60;
      kTraj = makeTrajPts( 2, 'propeller', nReadout, nLines, dkLine, nAngles );
      %nSpokes = 400;
      %nPtsPerSpoke = 100;
      %kTraj = makeTrajPts( 2, 'radial', nSpokes, nPtsPerSpoke );
      [fVals,phantImg] = makeNickPhantom( kTraj );
      N = [256, 256];
      psfMask = ones(2*N);

    case 1
      % Simulated data of NickPhantom with radial trajectory
      nSpokes = 400;
      nPtsPerSpoke = 100;
      kTraj = makeTrajPts( 2, 'radial', nSpokes, nPtsPerSpoke );
      [fVals,phantImg] = makeNickPhantom( kTraj );
      N = [256, 256];
      psfMask = ones(2*N);
      
    case 2
      % Simulated data of NickPhantom with rosette trajectory
      f1 = 147;
      f2 = 20;
      nK = 50000;
      dt = 0.001;
      kTraj = makeTrajPts( 2, 'rosette', nK, dt, f1, f2 );
      [fVals,phantImg] = makeNickPhantom( kTraj );
      N = [256, 256];
      psfMask = ones(2*N);
      
    case 3
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

    case 4
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

    case 5
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

  %figure; scatter( kTraj(:,1), kTraj(:,2), 4, 'k', 'filled' );  axis('equal');
end
