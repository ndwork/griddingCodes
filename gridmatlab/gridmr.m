function m = gridmr( d, k, n, varargin )
%function [m,nkx,nkxt,kwx] = gridmr( d, k, n, varargin ) % DEBUG
%GRIDMR  Puts k-space data on cartesian grid for inverse-FFT reconstruction
% Usages are
%   M = GRIDMR( D, K, N ) interpolates with triangle kernel
%   input:
%       d -- k-space data
%       k -- k-trajectory, scaled -0.5 to 0.5
%       n -- image size
%   output:
%       m -- gridded k-space data (n X n)
%   
%   M = GRIDMR( D, K, N, S ) interpolates with triangle kernel
%   input:
%       d -- k-space data
%       k -- k-trajectory, scaled -0.5 to 0.5
%       n -- image size
%       s -- gridding kernel size (number of grid points; W in Beatty '05)
%   output:
%       m -- gridded k-space data (n X n)
%
%   M = GRIDMR( D, K, N, S, A ) interpolates with kaiser-bessel kernel
%   computed with optimal shape parameter determined by A and S.
%   input:
%       d -- k-space data
%       k -- k-trajectory, scaled -0.5 to 0.5
%       n -- image size
%       s -- gridding kernel size (number of grid points; W in Beatty '05)
%       a -- grid oversampling ratio (determines kaiser-bessel shape)
%   output:
%       m -- gridded k-space data (n X n)
%
% Note: images from ifft of the gridded k-space data are apodised!  The de-
% apodisation function for a triangle kernel is
%    da(nx,ny) = (sinc(nx/N) .* sinc(ny/N)).^(-2)
% and the function for a kaiser-bessel kernel is
%    da(nx,ny) = (sinc(sqrt((s/(al*N)*nx).^2-(be/pi)^2))
%                 .* sinc(sqrt((s/(al*N)*ny).?2-(be/pi)^2).^(-1)
% (both over range nx,ny = -al*n/2, ..., al*n/2-1).  If the entire FOV of
% the a-oversampled image is de-apodised, regularisation may be needed to
% avoid divide-by-zero.
%
% Ethan Johnson, 2011

t = 4; % width of kernel in units of k-space samples (t = s-1)

% for alternate-width kernels
if( nargin >= 4 ); t = varargin{1}-1; end;

% for kaiser-bessel kernel
if( nargin >= 5 && ~isempty(varargin{2}) );
    al = varargin{2};
    be = (pi * sqrt( ((t+1)/al*(al-0.5))^2 - 0.8 ));
else
    be = -1;
end;


% convert to single column
d = d(:);
k = k(:);

% convert k-space samples to matrix-index units (not integers, though)
nkx = (n/2+1) + n*real(k);
nky = (n/2+1) + n*imag(k);

% allocate
m = zeros(n,n);

% loop over samples in kernel
for lx = -0.5*(t-1) : 1 : 0.5*(t-1)
% for lx = -t/2:t/2 % bad: redundant computations
    
    % find nearest samples
    nkxt = round(nkx+lx);
        
    % compute weighting from kernel
    if( be < 0 ); kwx = max(1-abs(nkx-nkxt)*(2/t),0);
    else kwx = kaibes( nkx-nkxt, t+1, be ); end;
    
    % map samples outside the matrix to the edges (suboptimal??)
    nkxt = max(nkxt,1); nkxt = min(nkxt,n);
    
    for ly = -0.5*(t-1) : 1 : 0.5*(t-1)
%     for ly = -t/2:t/2 % bad: redundant computations
        
        % find nearest samples
        nkyt = round(nky+ly);
        
        % compute weighting from kernel
        if( be < 0 ); kwy = max(1-abs(nky-nkyt)*(2/t),0);
        else kwy = kaibes( nky-nkyt, t+1, be ); end;
        
        % map samples outside the matrix to the edges
        nkyt = max(nkyt,1); nkyt = min(nkyt,n);
        
        % use sparse matrix to turn k-space trajectory into 2D matrix
        m = m+sparse(nkxt,nkyt,d.*kwx.*kwy,n,n);
        
    end
end

% zero out edge samples, since these may be due to samples outside
% the matrix
m(:,1) = 0; m(:,n) = 0;
m(1,:) = 0; m(n,:) = 0;

end

function c = kaibes( n, s, be )
%KAIBES

c = besseli( 0, be * sqrt( 1 - (2*n/s).^2 ) ); % in terms of s, beta
c( abs(n) > s/2 ) = 0; % should never be needed if gridding loop is correct

end