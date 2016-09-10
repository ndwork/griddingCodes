
function makePaperImgs_fp
  close all; clear;

  datacases = [2,1];
  downsamples = [4];
  outDir = './paperImgs_fp';


  mkdir( outDir );
  outFile = [outDir, '/makePaperImgsOut_fp.csv'];
  outFileID = fopen( outFile, 'w' );
  fprintf( outFileID, 'Datacase, Algorithm, Downsample, nIterations, MSE\n' );

  nIterations = [1 10 20 50 100 500 1000 2000 5000 10000];

  for datacaseIndx=1:numel(datacases)
    datacase = datacases(datacaseIndx);
    datacaseDir = [outDir, '/datacase_', num2str(datacase,'%4.4i')];
    mkdir( datacaseDir );

    for dsIndx=1:numel(downsamples)
      ds = downsamples(dsIndx);
      thisDs = num2str(ds, '%3.3i');

      [kTraj,iGridFVals,N] = loadDataCase( datacase, ds );
      nCoils = size( iGridFVals, 2 );

      for nIterIndx=1:numel(nIterations)
        nIter = nIterations(nIterIndx);
        thisNIter = num2str(nIter, '%4.4i');

        disp([datestr(now), ',  Working on fp with ', num2str(nIter), ...
          ' iterations and ds: ', thisDs ]);

        [weights,~,res] = makePrecompWeights_2D( kTraj, N, 'alg', 'fp', ...
          'nIter', nIter );
        close;

        if datacase==1
          cgrid = 2*N;
          radialImg = makeRadialImg( cgrid );
          mask = double( radialImg <= min(N) );
        else
          mask = [];
        end
        psfFile = [ datacaseDir, '/psf_', thisDs, '_fp_', thisNIter, '.jpg' ];
        [~,mse] = showPSF( weights, kTraj, N, mask, 'savefile', psfFile );

        figH = figure;
        scatter( kTraj(:,1), kTraj(:,2), 5, weights, 'filled' );
        colormap('copper');
        axis('equal');
        set( gca, 'xtick', [], 'ytick', [] );
        title([ 'Min=', num2str(min(weights)), ', Max=', num2str(max(weights)) ]);
        outFile = [ datacaseDir, '/paperImg_', thisDs, '_fp_', thisNIter, '_weights.png' ];
        saveas( figH, outFile );

        gridded = zeros( [ N nCoils] );
        for coilIndx=1:nCoils
          disp(['Working on coil:   ', num2str(coilIndx)]);
          gridded_2D = grid_2D( iGridFVals(:,coilIndx), kTraj, N, weights );
          gridded(:,:,coilIndx) = gridded_2D;
        end
        recon = sqrt( sum( gridded.*conj(gridded), 3 ) );

        scaleBounds = [ mean(recon(:)) - 2.5*std(recon(:)), ...
                      mean(recon(:)) + 2.5*std(recon(:))  ];
        scaledRecon = scaleImg( recon, scaleBounds );
        outFile = [ datacaseDir, '/paperImg_', thisDs, '_fp_', thisNIter, '.png' ];
        imwrite( scaledRecon, outFile );

        save( [datacaseDir,'/weights_', thisDs, '_fp_', thisNIter, '.mat'], ...
          'weights' );
        fprintf( outFileID, ['%4.4i fp, ', thisDs, thisNIter, ', %12.10f \n'], ...
          datacase, mse );

        close all;
        drawnow;
      end
    end
  end

  fclose( outFileID );
end

