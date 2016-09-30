
function makeFigures
  close all; clear;

  inDir = './paperImgs';
  outDir = './paperFigs';
  dataDir = '/Users/ndwork/Desktop/griddingData/';

  
  datacaseDirs = dir( [inDir,'/datacase_*'] );
  mkdir( outDir );
  
  for datacase = 1:numel(datacaseDirs)

    switch datacase
      case 1
        load([dataDir,'/data_2drad_test.mat']);
        ds = 4;
        N = [256 256];
      case 2
        load([dataDir,'/data_2dspiralNav.mat']);
        ds = 1;
        N = [128 64];
    end

    thisInDir = [ inDir, '/', datacaseDirs(datacase).name ];
    thisOutDir = [ outDir, '/', datacaseDirs(datacase).name ];
    mkdir( thisOutDir );


    % make weight images
    wFiles = dir( [thisInDir, '/weights_*.mat'] );
    ks = k(1:end,1:ds:end);
    ks = ks(:);
    kTraj = [ real(ks) imag(ks) ];
    nCoils = size(d,3);
    for i=1:numel(wFiles)
      name = wFiles(i).name;
      load( [thisInDir, '/', name] );
      figH = figure;
      scatter( kTraj(:,1), kTraj(:,2), 5, weights, 'filled' );
      colormap('copper');  caxis([-1d-4 1d-4]);
      axis('equal');
      set( gca, 'xtick', [], 'ytick', [] );
      nameParts = strsplit( name, {'_','.'} );
      if ~exist('nParts'), nParts = numel(nameParts); end;
      outFile = [ thisOutDir, '/', strjoin(nameParts(1:end-1),'_'), '.png' ];
      img = frame2im(getframe(figH));
      img = img(65:746,235:930,:);
      imwrite( img, outFile );

      reconNameParts = strsplit( name, {'_','.'});
      reconNameParts{1} = 'paperImg';
      if numel(reconNameParts) == nParts+1
        reconNameParts = { reconNameParts{1}, ...
          strjoin( reconNameParts(2:3), '_' ), reconNameParts{4:end} };
      end
      if numel( reconNameParts ) > 3
        tmp = reconNameParts{2};
        reconNameParts{2} = reconNameParts{3};
        reconNameParts{3} = tmp;
      end
      reconName = [ strjoin( reconNameParts(1:end-1), '_' ), '.png' ];
      reconImg = imread( [thisInDir, '/', reconName] );
      N = size(reconImg);
      gridded = zeros( [ N nCoils] );
      for coilIndx=1:nCoils
        iGridFVals = d(1:end,1:ds:end,coilIndx);
        iGridFVals = iGridFVals(:);
        gridded_2D = grid_2D( iGridFVals, kTraj, N, weights );
        gridded(:,:,coilIndx) = gridded_2D;
      end
      recon = sqrt( sum( gridded.*conj(gridded), 3 ) );
      scaledImg = scaleImg( recon, [0 0.3] );
      scaledNameParts = reconNameParts;
      scaledNameParts{1} = 'scaled';
      scaledName = [ strjoin( scaledNameParts(1:3), '_' ), '.png' ];
      imwrite( scaledImg, [thisOutDir, '/', scaledName] );

      close all; drawnow;
    end

    % Make noise sensitivity plot
    ks = k(1:end,1:ds:end);
    ks = ks(:);
    kTraj = [ real(ks) imag(ks) ];
    nTraj = size( kTraj, 1 );
    wFiles = dir( [thisInDir, '/weights_*.mat'] );
    for i=1:numel(wFiles)
      rlsdcMatch = regexp(wFiles(i).name, 'LSDC_con', 'once');
      lsdcMatch = regexp(wFiles(i).name, 'LSDC', 'once');
      if ~isempty( rlsdcMatch )
        rlsdcFile = wFiles(i).name;
      elseif ~isempty( lsdcMatch )
        lsdcFile = wFiles(i).name;
      end
    end
    load( [thisInDir, '/', lsdcFile] );
    lsdcWeights = weights;  clear( 'weights' );
    load( [thisInDir, '/', rlsdcFile] );
    rlsdcWeights = weights; clear( 'weights' );

    noiseScales = 0:0.2:1.0;
    nNoises = numel(noiseScales);
    lsdcMSEs = zeros( nNoises, 1 );
    rlsdcMSEs = zeros( nNoises, 1 );
    nWeights = numel(lsdcWeights);
    alpha = 1.5;
    W = 8;
    nC = 500;
    nGrid = ceil( alpha * N );
    trueAlpha = max( nGrid ./ N );
    [kCy,Cy,~] = makeKbKernel( N(1), N(1), trueAlpha, W, nC );
    [kCx,Cx,~] = makeKbKernel( N(2), N(2), trueAlpha, W, nC );
    for noiseIndx = 1:nNoises
      disp( ['Working on noiseIndx ', num2str(noiseIndx), ' of ', ...
        num2str(nNoises)] );
      noiseScale = noiseScales(noiseIndx);
      noise = randn(nWeights,1) * noiseScale;
      lsdcNoisyWeights = lsdcWeights .* (1+noise);
      rlsdcNoisyWeights = rlsdcWeights .* (1+noise);
      lsdcCw = applyC_2D( lsdcNoisyWeights, kTraj, 0, kCy, kCx, Cy, Cx, kTraj );
      lsdcCw = lsdcCw / median(lsdcCw);
      lsdcMSE = norm( lsdcCw(:) - 1, 2 ) / nWeights;
      rlsdcCw = applyC_2D( rlsdcNoisyWeights, kTraj, 0, kCy, kCx, Cy, Cx, kTraj );
      rlsdcCw = rlsdcCw / median(rlsdcCw);
      rlsdcMSE = norm( rlsdcCw(:) - 1, 2 ) / nWeights;
      lsdcMSEs(noiseIndx) = lsdcMSE;
      rlsdcMSEs(noiseIndx) = rlsdcMSE;
      close all;  drawnow;
    end
    noiseH = figure;
    plot( noiseScales, lsdcMSEs, 'k--', 'lineWidth', 2 );
    hold on;
    plot( noiseScales, rlsdcMSEs, 'b', 'lineWidth', 2 );
    set( gca, 'xtick', [0 1] );
    ax = gca;
    set( gca, 'ytick', [min(ax.YTick) max(ax.YTick)] );
    set(gca,'fontsize',20);
    legend( 'LSDC', 'rLSDC' );
    saveas( noiseH, [thisOutDir, '/noiseSensitivity.png'] );

  end
end

