
function [kTraj,iGridFVals,N] = loadDataCase( datacase, varargin )
  % [kTraj,iGridFVals] = loadDataCase( datacase [, ds ] )
  %
  % ds is the amount to downsample the data

  dataDir = '/Volumes/Seagate2TB/Data/griddingData/';
  %dataDir = '/Users/ndwork/Desktop/griddingData/';

  p = inputParser;
  p.addOptional( 'ds', 1, @isnumeric );
  p.parse(varargin{:});
  ds = p.Results.ds;

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

      N = [ 128 64 ];

  end

end