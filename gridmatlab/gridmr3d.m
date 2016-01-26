function m = gridmr3d( d, k, n, varargin )
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
%    da(nx,ny) = (sinc(nx/N) .* sinc(ny/N) .* sinc(nz/N)).^(-2)
% and the function for a kaiser-bessel kernel is
%    da(nx,ny) = (sinc(sqrt((s/(al*N)*nx).^2-(be/pi)^2))
%                 .* sinc(sqrt((s/(al*N)*ny).^2-(be/pi)^2).^(-1)
% (both over range nx,ny = -al*n/2, ..., al*n/2-1).  If the entire FOV of
% the a-oversampled image is de-apodised, regularisation may be needed to
% avoid divide-by-zero.
%
% Ethan Johnson, 2011

t = 4; % width of kernel in units of k-space samples (t = s-1)

% for alternate-width kernels
if( nargin >= 4 && ~isempty(varargin{1}) ); t = varargin{1}-1; end;

% for kaiser-bessel kernel
if( nargin >= 5 && ~isempty(varargin{2}) );
    al = varargin{2};
    be = (pi * sqrt( ((t+1)/al*(al-0.5))^2 - 0.8 ));
    fprintf('*** kaiser-bessel. ***\n'); % DEBUG
else
    be = -1;
    fprintf('*** triangle. ***\n'); % DEBUG
end;

% for variable-resolution/field-of-view acquisitions
if isscalar(n), nx = n; ny = n; nz = n;
elseif length(n)==3, nx = n(1); ny = n(2); nz = n(3);
else nx = n(1); ny = n(1); nz = n(1); % fallback position
end
 

% convert to single column
d = d(:);
k = reshape(k,[],3);

% convert k-space samples to matrix-index units (not integers, though)
nkx = (nx/2+1) + nx*k(:,1);
nky = (ny/2+1) + ny*k(:,2);
nkz = (nz/2+1) + nz*k(:,3);

% allocate
m = zeros(nx,ny*nz); % reshape to 3d later

% loop over samples in kernel
for lx = -0.5*(t-1) : 1 : 0.5*(t-1)
% for lx = -t/2:t/2 % bad: redundant computations
    
    % find nearest samples
    nkxt = round(nkx+lx);
        
    % compute weighting from kernel
    if( be < 0 ); kwx = max(1-abs(nkx-nkxt)*(2/t),0); % triangle
    else kwx = kaibes( nkx-nkxt, t+1, be ); end; % kaiser-bessel
    
    % map samples outside the matrix to the edges (suboptimal??)
    nkxt = max(nkxt,1); nkxt = min(nkxt,nx);
    
    for ly = -0.5*(t-1) : 1 : 0.5*(t-1)
%     for ly = -t/2:t/2 % bad: redundant computations
        
        % find nearest samples
        nkyt = round(nky+ly);
        
        % compute weighting from kernel
        if( be < 0 ); kwy = max(1-abs(nky-nkyt)*(2/t),0); % triangle
        else kwy = kaibes( nky-nkyt, t+1, be ); end; % kaiser-bessel
        
        % map samples outside the matrix to the edges
        nkyt = max(nkyt,1); nkyt = min(nkyt,ny);
        
        for lz = -0.5*(t-1) : 1 : 0.5*(t-1)
            
            % find nearest samples
            nkzt = round(nkz+lz);
            
            % compute weighting from kernel
            if( be < 0 ); kwz = max(1-abs(nkz-nkzt)*(2/t),0); % triangle
            else kwz = kaibes( nkz-nkzt, t+1, be ); end; % kaiser-bessel
            
            % map samples outside the matrix to the edges
            nkzt = max(nkzt,1); nkzt = min(nkzt,nz);

            % ----- (the punch line!):
            if( sum( abs((nkyt+(nkzt-1)*ny)-round(nkyt+(nkzt-1)*ny))>0 )>0 ),
              fprintf('*** problem?!\n');  % DEBUG
            end
            m = m+sparse(nkxt,nkyt+(nkzt-1)*ny,d.*kwx.*kwy.*kwz,nx,ny*nz);
            % -----
        end
    end
end

% reshape to 3d
m = reshape(m,[nx ny nz]);

% zero out edge-face samples, since these may be due to samples outside
% the matrix
m(:,:,1) = 0; m(:,:,nz) = 0;
m(:,1,:) = 0; m(:,ny,:) = 0;
m(1,:,:) = 0; m(nx,:,:) = 0;

end

function c = kaibes( n, s, be )
%KAIBES

c = besseli( 0, be * sqrt( 1 - (2*n/s).^2 ) ); % in terms of s, beta
c( abs(n) > s/2 ) = 0; % should never be needed if gridding loop is correct

end