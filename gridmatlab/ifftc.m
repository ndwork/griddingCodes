function g = ifftc( f, varargin )
%IFFTC 'centred' Fourier transform
%   IFFTC( X ) returns the fourier transform of X with the spectrum centred
%   in the output, presuming X is centered in the input.  The 'symmetric',
%   'nonsymmetric' stuff is not supported.
%   Essentially it is IFFT with FFTSHIFTs.  'Centred' means that the (0,0)
%   coordinate is at X(floor(end/2)+1).

if( nargin == 3 ); dim = varargin{2};
% elseif( size(f,1) == 1 && size(f,2) > 1 ); dim = 2;
% else dim = 1;
elseif numel(f)==1, dim = 1;
else dim = find(size(f)>1,1);
end

if( nargin == 1 )
    g = fftshift( ifft( ifftshift( f, dim ), [], dim ), dim );
elseif( nargin < 4 )
    g = fftshift( ifft( ifftshift( f, dim ), varargin{1}, dim ), dim );
else
    error( 'Too many input arguments' );
end

end