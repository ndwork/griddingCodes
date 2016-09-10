
function k = iGrid( img, traj, varargin )
  % k = applyE( img, traj, [ 'alpha', alpha, 'w', w, 'nC', nC ] )
  % MRI encoding with Inverse Gridding
  % img is a 2D array specifying the volume to be encoded
  % traj is a MxV array specifying the k-space trajectory.
  %   V is the number of dimensions of the img
  %   The first/second/third column is kx/ky/kz
  %   The units are normalized to [-0.5,0.5).
  % k is a 1D array of M elements specifying the k-space data values
    % alpha is the oversampling factor [1,inf]
  % W is the window width in pixels
  % nC is the number of points to sample the convolution kernel
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on Beatty et. al., IEEE TMI, 2005

  if size( traj, 1 ) == 1, traj=transpose(traj); end;

  nD = size( traj, 2 );
  if nD == 1
    out = iGrid_1D( img, traj, varargin{:} );
  elseif nD == 2
    k = iGrid_2D( img, traj, varargin{:} );
  elseif
    k = iGrid_3D( img, traj, varargin{:} );
  end
end
