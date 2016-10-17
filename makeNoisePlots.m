
function makeNoisePlots
  close all; clear;

  datacases = [1,2];
  outDir = './paperImgs';
  algorithms = {'LSDC','rLSDC'};
  Ns = { [ 256 256 ], [ 128 64 ] };

  noiseStdDevs = 0:0.1:1;

  lineStyles = {'-','--',':'};
  mses = zeros(numel(algorithms),numel(noiseStdDevs));
  for datacaseIndx=1:numel(datacases)
    datacase = datacases(datacaseIndx);
    datacaseDir = [outDir, '/datacase_', num2str(datacase,'%4.4i')];

    N = Ns{datacaseIndx};

    ellipticalImg = makeEllipticalImage( 2*N );
    mask = double( ellipticalImg <= 0 );

    for noiseIndx = 1:numel(noiseStdDevs)

      for algIndx=1:numel(algorithms)
        thisAlg = algorithms{algIndx};

        load( [datacaseDir,'/weights_',thisAlg,'.mat'] );
        nTraj = size(kTraj,1);

        noise = noiseStdDevs(noiseIndx) * randn( nTraj, 1 );
        noisyWeights = weights .* ( 1 + noise );
        [~,mse] = showPSF( noisyWeights, kTraj, N, mask ); close;
        mses( algIndx, noiseIndx ) = mse;
        drawnow;
      end
    end

    figure; hold all;
    for algIndx=1:numel(algorithms)
      thisLineStyle = lineStyles{mod(algIndx-1,numel(lineStyles))+1};
      plot( noiseStdDevs, mses(algIndx,:), 'LineWidth', 3, ...
        'LineStyle', thisLineStyle);
    end
    legendH = legend( algorithms{:} );
    set( legendH, 'FontSize', 20, 'Position', [0.2 0.75 0.2 0.1] );
    ax = gca;
    ax.FontSize = 20;
    ax.XTick = [min(noiseStdDevs) max(noiseStdDevs)];
    ax.YTick = [min(mses(:)) max(mses(:))];
    ax.YAxis.TickLabelFormat = '%,.1f';
    ax.LineWidth = 1.5;
    drawnow;

  end

end
