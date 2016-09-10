
function makePaperSims
  close all;  clear;

  % Function Parameters
  nPts = 64;
  trajs = {'random', 'poissonDisc', 'radial'};
  algs = {'sdLSDC','cls','con_sdLSDC','fdLSDC','fp'};
  outFile = 'makePaperSimsOut.csv';
  outDir = 'paperSims';
  noiseFraction = 0.3;


  % Code
  N = [nPts nPts];
  trajParams = containers.Map();
  trajParams('radial') = {360,50};
  trajParams('random') = {round(nPts*nPts)};
  trajParams('poissonDisc') = {1/nPts};

  mkdir( outDir );
  outFileID = fopen( [outDir,'/',outFile], 'w' );
  fprintf( outFileID, 'Trajectory, Algorithm, MSE, noisyMSE, min weight, max weight\n' );
  fclose( outFileID );


  for trajIndx = 1:numel(trajs)
    thisTraj = trajs{trajIndx};
    thisTrajParam = trajParams(trajs{trajIndx});
    kTraj = makeTrajPts( 2, thisTraj, thisTrajParam{:} );

    trajFigH = figure;
    scatter( kTraj(:,1), kTraj(:,2), 8, 'ko', 'fill' );
    set( gca, 'xtick', [], 'ytick', [] );
    axis('equal');
    trajFile = [outDir, '/', thisTraj,'.jpg'];
    saveas( trajFigH, trajFile );

    for algIndx = 1:numel(algs)
      thisAlg = algs{algIndx};
      weights = makePrecompWeights_2D( kTraj, N, 'alg', thisAlg );

      cgrid = 2*N;
      radialImg = makeRadialImg( cgrid );
      mask = double( radialImg <= min(N) );
      psfFile = [ outDir, '/', thisTraj, '_', thisAlg, '_psf.jpg' ];
      [~,mse] = showPSF( weights, kTraj, N, mask, 'savefile', psfFile );

      figH = figure;
      scatter( kTraj(:,1), kTraj(:,2), 5, weights, 'filled' );
      colormap('copper');
      axis('equal');
      set( gca, 'xtick', [], 'ytick', [] );
      title([ 'Min=', num2str(min(weights)), ', Max=', num2str(max(weights)) ]);
      weightFile = [ outDir, '/', thisTraj, '_', thisAlg, '_weights.jpg' ];
      saveas( figH, weightFile );

      noise = randn(size(weights))*noiseFraction;
      noisyWeights = weights .* ( 1 + noise );
      psfFile = [ outDir, '/', thisTraj, '_', thisAlg, '_psf_noisy.jpg' ];
      [~,noisyMSE]=showPSF( noisyWeights, kTraj, N, mask, ...
        'savefile', psfFile );

      outFileID = fopen( [outDir,'/',outFile], 'a' );
      fprintf( outFileID, [thisTraj, ', ', thisAlg, ...
        ', %14.10f, %14.10f, %14.10f, %14.10f \n'], mse, noisyMSE, ...
        min(weights(:)), max(weights(:)) );
      fclose( outFileID );

      close all;
    end
  end

end
