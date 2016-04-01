
function out = iGridT_2D( F, traj, N, varargin )
  % out = iGridT_2D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
  %   traj is a Mx2 array specifying the k-space trajectory.
  %     The first/second column is kx/ky
  %     The units are normalized to [-0.5,0.5).
  %   N is the size of the output image
  %     If N is a scalar, then the final image is assumed to be square
  %     If N is a 2 element array, then N = [Ny Nx]
  %
  % Optional Inputs:
  %   alpha - a float parameter specifying the oversampling factor
  %   W - an integer specifying the kernel's width
  %   nC - specifies the number of samples in the kernel
  %   verbose - if set to true, outputs processing status
  %
  % Written by Nicholas Dwork (c) 2015
  % Based on EE369C notes written by John Pauly

  if numel(N)==1
    Ny=N;  Nx=N; 
  else
    Ny=N(1); Nx=N(2);
  end

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParamValue( 'alpha', defaultAlpha, checknum );
  p.addParamValue( 'W', defaultW, checknum );
  p.addParamValue( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  nGridY = Ny;
  nGridX = Nx;

  %% Make the convolution kernel
  Gy = nGridY;
  [kCy,Cy,cImgY,kwy] = makeKbKernel( Gy, Ny, alpha, W, nC );
  Gx = nGridX;
  [kCx,Cx,cImgX,kwx] = makeKbKernel( Gx, Nx, alpha, W, nC );
  kws = [ kwy kwx ];

  gridKs = size2fftCoordinates([ nGridY nGridX ]);
  gridKy=gridKs{1};  gridKx=gridKs{2};
  [gridKx,gridKy] = meshgrid(gridKx,gridKy);

  % Perform Adjoint of Circular Convolution
  nTraj = size(traj,1);
  kDistThreshY = 0.5*kwy;
  kDistThreshX = 0.5*kwx;
  fftGridded = zeros( nGridY, nGridX );
  for trajIndx=1:nTraj
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThreshY & distsKx < kDistThreshX );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    fftGridded(shortDistIndxs) = fftGridded(shortDistIndxs) + ...
      F(trajIndx) * ( CValsY .* CValsX );
  end

  % Circular convolution
  onesCol = ones(nTraj,1);
  for dim=1:2
    alt = zeros( size(traj) );
    alt(:,dim) = onesCol;

    for altDir=[-1 1]
      NewTraj = traj + altDir*alt;
      if altDir < 0
        NewTrajIndxs = find( NewTraj(:,dim) > -0.5-kws(dim)/2 );
      else
        NewTrajIndxs = find( NewTraj(:,dim) < 0.5+kws(dim)/2 );
      end

      NewTraj = NewTraj( NewTrajIndxs, : );
      for i=1:numel(NewTrajIndxs)
        trajIndx = NewTrajIndxs(i);
        NewDistsKy = abs( NewTraj(i,1) - gridKy );
        NewDistsKx = abs( NewTraj(i,2) - gridKx );
        NewShortDistIndxs = find( NewDistsKy < kDistThreshY & ...
                                  NewDistsKx < kDistThreshY );
        NewShortDistsKy = NewDistsKy( NewShortDistIndxs );
        NewShortDistsKx = NewDistsKx( NewShortDistIndxs );
        NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
        NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
        fftGridded(NewShortDistIndxs) = fftGridded(NewShortDistIndxs) + ...
          F(trajIndx) * ( NewCValsY .* NewCValsX );
      end
    end
  end

  % Perform an inverse fft
  data = fftshift( ifft2( ifftshift(fftGridded) ) );

  % Perform deapodization
  cImg = transpose(cImgY) * cImgX;
  out = data ./ cImg;
end

