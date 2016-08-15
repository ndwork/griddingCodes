
function [scale,mse] = showPSF( weights, traj, N, varargin )
  % scale = showPSF( weights, traj, N, mask [, ...
  %   'savefile', savefile, 'imgTitle', imgTitle ] )
  %
  % Written by Nicholas Dwork - Copyright 2016

  p = inputParser;
  p.addOptional( 'mask', [] );
  p.addParamValue( 'savefile', [] );
  p.addParamValue( 'imgTitle', [] );
  p.parse( varargin{:} );
  mask = p.Results.mask;
  savefile = p.Results.savefile;
  imgTitle = p.Results.imgTitle;

  nGrid = 2 * N;
  W = 8;
  nC = 500;

  psf = iGridT_2D( weights, traj, nGrid, 'alpha', 2.0, 'W', W, 'nC', nC );
  psf = cropData( psf, 2*N );

  if isempty(mask)
    %radialImg = makeRadialImg( size(psf) );
    %mask = radialImg <= min(N);
    mask = ones( size(psf) );
  end
  
  mask = cropData( mask, 2*N );
  psf = mask .* psf;
  absPsf=abs(psf);
  absPsf_dB = 20*log10(absPsf);
  absPsf_dB(~isfinite(absPsf_dB)) = 0;

  scale = 1 ./ max(absPsf(:));
  
  figure;
  imshow( imresize(absPsf_dB,2,'nearest'), [-100 0] );
  if ~isempty(imgTitle), title(imgTitle); end;
  drawnow;

  if ~isempty(savefile)
    scaled_dB = scaleImg( absPsf_dB, [-100 0], [0 1] );
    imwrite( scaled_dB, savefile );
  end
  
  b=zeros(size(psf)); b(1,1)=1; b=fftshift(b);
  mse = sum( abs( mask(:).*psf(:) - mask(:).*b(:) ).^2 ) ./ sum(mask(:));
  disp(['MSE: ', num2str(mse)]);
end
