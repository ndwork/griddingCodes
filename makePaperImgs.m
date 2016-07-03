
function makePaperImgs

  outDir = './paperImgs';

  dataDir = '/Volumes/Seagate2TB/Data/griddingData/';
  dataFile = 'data_2drad_test.mat';

  load( [dataDir dataFile] );

  % k contains the trajectory data
  % each column is a radial line

  downsamples = [ 16 8 4 2 ];
  algorithms = { 'cls', 'fp', 'ls', 'rls' };

  mkdir( outDir );

  N = [ 256 256 ];
  for dsIndx=1:numel(downsamples)
    ds = downsamples(dsIndx);

    iGridFVals = d(:,1:ds:end);
    iGridFVals = iGridFVals(:);

    ks = k(:,1:ds:end);
    ks = ks(:);
    kTraj = [ real(ks) imag(ks) ];

    for algIndx=1:numel(algorithms)
      thisAlg = algorithms{algIndx};
      thisDs = num2str(ds, '%3.3i');

      disp([ 'Working on ', thisAlg, ' algorithm with ds: ', thisDs ]);

      weights = makePrecompWeights_2D( kTraj, N, 'alg', thisAlg );
      weights = weights / sum( weights(:) );
      gridded_2D = grid_2D( iGridFVals, kTraj, N, weights );

      recon = abs( gridded_2D );
      scaledRecon = scaleImg( recon );

      outFile = [ outDir, '/paperImg_', thisDs, '_', thisAlg, '.png' ];
      imwrite( scaledRecon, outFile );
      
      save( [outDir,'/weights_',thisAlg,'_',thisDs,'.mat'], 'weights' );
    end
  end

end

