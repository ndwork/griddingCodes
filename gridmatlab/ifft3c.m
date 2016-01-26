function g = ifft3c( f, varargin )
% IFFT3C Three-dimentional inverse 'centred' Fourier transform
%   IFFT3C( X ) computes the three-dimensional inverse fourier transform of
%   X with the spectrum centred in the output, presuming that X is centered
%   in the input (essentially it is IFFT3 with FFTSHIFTs).  'Centred' means
%   that the (0,0) coordinate is at index floor(end/2)+1 in the first three
%   dimensions.  (There is no actual 'IFFT3' function.)
%
% Ethan Johnson, 2013

sz = size(f);

if nargin>1 && length(varargin{1})==length(sz), sz = varargin{1}; end

g = ifftc(ifftc(ifftc( f, sz(1),1),sz(2),2),sz(3),3);

end