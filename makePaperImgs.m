
function makePaperImgs
  close all; clear;

  outFile = './paperImgs/makePaperImgsOut.csv';
  outDir = './paperImgs';

  %dataDir = '/Volumes/Seagate2TB/Data/griddingData/';
  dataDir = '/Users/ndwork/Desktop/';
  dataFile = 'data_2drad_test.mat';

  mkdir( outDir );
  outFileID = fopen( outFile, 'w' );
  fprintf( outFileID, 'Algorithm, Downsample, MSE\n' );

  load( [dataDir dataFile] );
  % k contains the trajectory data
  % each column is a radial line

  downsamples = [ 16 8 4 2 1 ];
  %algorithms = { 'cls', 'fp', 'regWLS' };
  algorithms = {'wls','fp'};

  N = [ 256 256 ];
  for dsIndx=1:numel(downsamples)
    ds = downsamples(dsIndx);

    iGridFVals = d(1:2:end,1:ds:end);
    iGridFVals = iGridFVals(:);

    ks = k(1:2:end,1:ds:end);
    ks = ks(:);
    kTraj = [ real(ks) imag(ks) ];

    for algIndx=1:numel(algorithms)
      thisAlg = algorithms{algIndx};
      thisDs = num2str(ds, '%3.3i');

      disp([datestr(now), ',  Working on ', thisAlg, ' algorithm with ds: ', thisDs ]);

      [weights,~,res] = makePrecompWeights_2D_v2( kTraj, N, 'alg', thisAlg );
      psfFile = [ outDir, '/psf_', thisDs, '_', thisAlg, '.jpg' ];
      [scale,mse] = showPSF( weights, kTraj, N, 'savefile', psfFile );
      weights = weights * scale;

      figH = figure;
      scatter( kTraj(:,1), kTraj(:,2), 5, weights, 'filled' );
      colormap('copper');
      axis('equal');
      title([ 'Min=', num2str(min(weights)), ', Max=', num2str(max(weights)) ]);
      outFile = [ outDir, '/paperImg_', thisDs, '_', thisAlg, '_weights.png' ];
      saveas( figH, outFile );

      gridded_2D = grid_2D( iGridFVals, kTraj, N, weights );
      recon = abs( gridded_2D );

      scaledRecon = scaleImg( recon );

      outFile = [ outDir, '/paperImg_', thisDs, '_', thisAlg, '.png' ];
      imwrite( scaledRecon, outFile );

      save( [outDir,'/weights_',thisAlg,'_',thisDs,'.mat'], 'weights' );
      fprintf( outFileID, [thisAlg, ', ', thisDs, ', %12.10f \n'], mse );

      close all;
      drawnow;
    end
  end

  fclose(outFileID);
end

