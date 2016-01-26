function g = fft3c( f, varargin )
% FFT3C Three-dimentional inverse 'centred' Fourier transform
%   FFT3C( X ) computes the three-dimensional Fourier transform of X with the spectrum centred in the output, presuming X is centered in the
%   input.  Essentially it is 'IFFT3' with FFTSHIFTs.  'Centred' means that
%   the (0,0) coordinate is at X(floor(end/2)+1,floor(end/2)+1).  (There is
%   no actual 'IFFT3' function)
%
% Ethan Johnson, 2013

sz = size(f);

if nargin>1 && length(varargin{1})==length(sz), sz = varargin{1}; end

g = fftc(fftc(fftc( f, sz(1),1),sz(2),2),sz(3),3);

end