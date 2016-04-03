
function out = iGridT_3D( F, traj, N, varargin )
  % out = iGridT_3D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
  %   traj is a Mx3 array specifying the k-space trajectory.
  %     The first/second/third column is kx/ky/kz
  %     The units are normalized to [-0.5,0.5).
  %   N is the size of the output image
  %     If N is a scalar, then the final image is assumed to be square
  %     If N is a 2 element array, then N = [Ny Nx Nz]
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
    Ny=N;  Nx=N;  Nz=N;
  else
    Ny=N(1); Nx=N(2); Nz=N(3);
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
  nGridZ = Nz;

  %% Make the convolution kernel
  Gy = nGridY;
  [kCy,Cy,cImgY,kwy] = makeKbKernel( Gy, Ny, alpha, W, nC );
  Gx = nGridX;
  [kCx,Cx,cImgX,kwx] = makeKbKernel( Gx, Nx, alpha, W, nC );
  Gz = nGridZ;
  [kCz,Cz,cImgZ,kwz] = makeKbKernel( Gz, Nz, alpha, W, nC );
  kws = [ kwy kwx kwz ];

  gridKs = size2fftCoordinates([ nGridY nGridX nGridZ ]);
  gridKy=gridKs{1};  gridKx=gridKs{2};  gridKz=gridKs{3};

  % Perform Adjoint of Circular Convolution
  nTraj = size(traj,1);
  kDistThreshY = 0.5*kwy;
  kDistThreshX = 0.5*kwx;
  kDistThreshZ = 0.5*kwz;
  fftGridded = zeros( nGridY, nGridX, nGridZ );
  for trajIndx = 1:nTraj
    distsKy = abs( traj(trajIndx,1) - gridKy );
    distsKx = abs( traj(trajIndx,2) - gridKx );
    distsKz = abs( traj(trajIndx,3) - gridKz );
    shortDistIndxsY = find( distsKy < kDistThreshY );
    shortDistIndxsX = find( distsKx < kDistThreshX );
    shortDistIndxsZ = find( distsKz < kDistThreshZ );
    shortDistsKy = distsKy( shortDistIndxsY );
    shortDistsKx = distsKx( shortDistIndxsX );
    shortDistsKz = distsKz( shortDistIndxsZ );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    CValsZ = interp1( kCz, Cz, shortDistsKz, 'linear', 0 );
    CVals = bsxfun( @times, CValsY*transpose(CValsX), ...
      reshape( CValsZ, [1 1 numel(CValsZ)] ) );
    fftGridded(shortDistIndxsY,shortDistIndxsX,shortDistIndxsZ ) = ...
      fftGridded(shortDistIndxsY,shortDistIndxsX,shortDistIndxsZ ) + ...
      F(trajIndx) * CVals;
  end

  % Circular convolution
  onesCol = ones(nTraj,1);
  for dim=1:3
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
        NewDistsKz = abs( NewTraj(i,3) - gridKz );
        NewShortDistIndxsY = find( NewDistsKy < kDistThreshY );
        NewShortDistIndxsX = find( NewDistsKx < kDistThreshX );
        NewShortDistIndxsZ = find( NewDistsKz < kDistThreshZ );
        NewShortDistsKy = NewDistsKy( NewShortDistIndxsY );
        NewShortDistsKx = NewDistsKx( NewShortDistIndxsX );
        NewShortDistsKz = NewDistsKz( NewShortDistIndxsZ );
        NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
        NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
        NewCValsZ = interp1( kCz, Cz, NewShortDistsKz, 'linear', 0 );
        NewCVals = bsxfun( @times, NewCValsY*transpose(NewCValsX), ...
          reshape( NewCValsZ, [1 1 numel(NewCValsZ)] ) );
        fftGridded( NewShortDistIndxsY, NewShortDistIndxsX, NewShortDistIndxsZ ) = ...
          fftGridded( NewShortDistIndxsY, NewShortDistIndxsX, NewShortDistIndxsZ ) + ...
          F(trajIndx) * NewCVals;
      end
    end
  end

  % Perform an inverse fft
  data = fftshift( ifftn( ifftshift(fftGridded) ) );

  % Perform deapodization
  cImgYX = transpose(cImgY) * cImgX;
  cImgZ_reshaped = reshape( cImgZ, [1 1 numel(cImgZ)] );
  cImg = bsxfun( @times, cImgZ_reshaped, cImgYX );
  out = data ./ cImg;
end

