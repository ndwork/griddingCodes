function g = fft2c( f, varargin )
%FFT2C Two-dimentional 'centred' Fourier transform
%   FFT2C( X ) returns the two-dimensional fourier transform of X with the
%   spectrum centred in the output, presuming X is centered in the input.
%   Essentially it is FFT2 with FFTSHIFTs.  'Centred' means that the (0,0)
%   coordinate is at X(floor(end/2)+1,floor(end/2)+1).

%g = fftshift( fft2( fftshift( f ), varargin{:} ) );
g = fftshift(fftshift(  fft2(  ifftshift(ifftshift( ...
    f ...
    ,1),2),  varargin{:})  ,1),2);

end