function g = fftc( f, varargin )
%FFTC 'centred' Fourier transform
%   FFTC( X ) returns the fourier transform of X with the spectrum centred
%   in the output, presuming X is centered in the input.
%   Essentially it is FFT with FFTSHIFTs.  'Centred' means that the (0,0)
%   coordinate is at X(floor(end/2)+1).
%   
%   Note: fftc( f, n ) doesn't do the right thing for n > size(f,1)
%   --fixed now, I think

if( nargin == 3 ); dim = varargin{2};
% elseif( size(f,1) == 1 && size(f,2) > 1 ); dim = 2;
% else dim = 1;
elseif numel(f)==1, dim = 1;
else dim = find(size(f)>1,1);
end

if( nargin == 1 )
    g = fftshift( fft( ifftshift( f, dim ), [], dim ), dim );
elseif( nargin < 4 )
    n = varargin{1};
    if isempty(n) || n < size(f,dim)
        g = fftshift( fft( ifftshift( f, dim ), [], dim ), dim );
%     elseif( ~isempty(n) && n < size(f,dim) ) %--2011/12/10: obviated by || modification above!
%         g = fftshift( fft( ifftshift( f, dim ), n, dim ), dim );
    else
        presz = zeros(1,ndims(f)); postsz = zeros(1,ndims(f));
        presz(dim) = ceil( (n - size(f,dim)) / 2 );
        postsz(dim) = floor( (n - size(f,dim)) / 2 );
        f = padarray( f, presz, 'pre' ); f = padarray( f, postsz, 'post' );
        g = fftshift( fft( ifftshift( f, dim ), [], dim ), dim );
    end
% elseif( nargin < 4 )
%     g = fftshift( fft( ifftshift( f, dim ), varargin{1}, dim ), dim );
else
    error( 'Too many input arguments' );
end

end