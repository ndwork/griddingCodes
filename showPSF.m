
function scale = showPSF( weights, traj, N, mask, varargin )
  % scale = showPSF( weights, traj, N, mask [, savefile ] )
  %
  % Written by Nicholas Dwork - Copyright 2016

  p = inputParser;
  p.addOptional( 'savefile', [] );
  p.parse( varargin{:} );
  savefile = p.Results.savefile;

  nGrid = 2.5 * N;
  W = 8;
  nC = 500;

  psf = iGridT_2D( weights, traj, nGrid, 'alpha', 2.0, 'W', W, 'nC', nC );
  psf = cropData( psf, 2*N );
  
  if nargin < 4
    radialImg = makeRadialImg( size(psf) );
    mask = radialImg <= min(N);
  end
  
  mask = cropData( mask, 2*N );
  psf = mask .* psf;
  psf = psf ./ max(psf(:));
  b=zeros(size(psf)); b(1,1)=1; b=fftshift(b);
  mse = sum( abs( mask(:).*psf(:) - mask(:).*b(:) ).^2 ) ./ sum(mask(:));
  disp(['MSE: ', num2str(mse)]);
  tmp=abs(psf);
  scale = 1 ./ max(tmp(:));
  tmp = tmp * scale;
  tmp=20*log10(tmp);
  tmp(~isfinite(tmp))=0;

  if ~isempty(savefile)
    F=getframe(gca);
    im=F.cdata;
    imwrite( savefile, im );
  end
  
  figure;
  imshow( imresize(tmp,2,'nearest'), [-100 0] );
  drawnow;
end
