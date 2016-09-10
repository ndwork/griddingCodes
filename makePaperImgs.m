
function makePaperImgs
  close all; clear;

  datacases = [2,1];
  ds = 4;
  outDir = './paperImgs';
  algorithms = {'fp','fdLSDC','fdLSDC_con','cls'};


  mkdir( outDir );
  outFile = [outDir, '/makePaperImgsOut.csv'];
  outFileID = fopen( outFile, 'a' );
  fprintf( outFileID, 'Datacase, Algorithm, MSE\n' );

  for datacaseIndx=1:numel(datacases)
    datacase = datacases(datacaseIndx);
    datacaseDir = [outDir, '/datacase_', num2str(datacase,'%4.4i')];
    mkdir( datacaseDir );

    [kTraj,iGridFVals,N] = loadDataCase( datacase, ds );
    nCoils = size( iGridFVals, 2 );
    
    for algIndx=1:numel(algorithms)
      thisAlg = algorithms{algIndx};

      disp([datestr(now), ',  Working on ', thisAlg ]);

      [weights,~,~] = makePrecompWeights_2D( kTraj, N, 'alg', thisAlg );
      close; drawnow;

      psfFile = [ datacaseDir, '/psf_', thisAlg, '.jpg' ];
      [~,mse] = showPSF( weights, kTraj, N, 'savefile', psfFile );
      save( [datacaseDir,'/weights_',thisAlg,'.mat'], 'weights', 'kTraj' );

      figH = figure;
      scatter( kTraj(:,1), kTraj(:,2), 5, weights, 'filled' );
      colormap('copper');
      axis('equal');
      set( gca, 'xtick', [], 'ytick', [] );
      title([ 'Min=', num2str(min(weights)), ', Max=', num2str(max(weights)) ]);
      outFile = [ datacaseDir, '/paperImg_', thisAlg, '_weights.png' ];
      saveas( figH, outFile );

      gridded = zeros( [ N nCoils] );
      for coilIndx=1:nCoils
        disp(['Working on traj/alg/coil:   ', thisAlg, ...
          ' / ', num2str(coilIndx)]);
        gridded_2D = grid_2D( iGridFVals(:,coilIndx), kTraj, N, weights );
        gridded(:,:,coilIndx) = gridded_2D;
      end
      recon = sqrt( sum( gridded.*conj(gridded), 3 ) );

      save( [datacaseDir, '/recon_', thisAlg, '.mat' ], 'recon' );

      scaleBounds = [ mean(recon(:)) - 2.5*std(recon(:)), ...
                      mean(recon(:)) + 2.5*std(recon(:))  ];
      scaledRecon = scaleImg( recon, scaleBounds );
      outFile = [ datacaseDir, '/paperImg_', thisAlg, '.png' ];
      imwrite( scaledRecon, outFile );

      fprintf( outFileID, ['%4.4i', thisAlg, ', %12.10f \n'], ...
        datacase, mse );

      close all;
      drawnow;
    end
  end

  fclose( outFileID );
end

