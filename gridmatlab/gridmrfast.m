function m = gridmrfast( d, k, n, varargin )
%GRIDMRFAST  Fast (for Kaiser-Bessel kernel) version of GRIDMR
% Usages are same as GRIDMR.
% Maximum error due to kernel sampling/interpolation is fixed below 1e-3.
%
% Ethan Johnson, 2011

t = 4; % (default) width of kernel in units of k-space samples (t = s-1)
if( nargin >= 4 ); t = varargin{1}-1; end; % for alternate-width kernel

% for kaiser-bessel kernel
if( nargin >= 5 );
    al = varargin{2};
    be = (pi * sqrt( ((t+1)/al*(al-0.5))^2 - 0.8 ));
    
    % --- nearest neighbour interpolation ---
%     l = ceil(0.91/(al*( 1e-3 ))); % kernel sampling density
%     nkl = ( -(t+1)/2 : 1/l : (t+1)/2-(1/l) )/(t+1); % kernel sample points
%     c = [ 0 besseli( 0, be*sqrt(1-(2*nkl).^2 )) 0 ].';
    
    % --- bilinear interpolation ---
    l = ceil(sqrt(0.37/(al*( 1e-3 )))); % kernel sampling density
    nkl = ( -(t+1)/2 : 1/l : (t+1)/2-(1/l) )/(t+1); % kernel sample points
    c = [ 0 0 besseli( 0, be*sqrt(1-(2*nkl).^2) ) 0 0 ].';

else
    l = -1; c = [];
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
    
    % find nearest samples
    nkxt = round(nkx+lx);
        
    % compute weighting from kernel
    if( l < 0 ); kwx = max(1-abs(nkx-nkxt)*(2/t),0);
    else kwx = kaibes( nkx-nkxt, t+1, c, l ); end;
    
    % map samples outside the matrix to the edges (suboptimal??)
    nkxt = max(nkxt,1); nkxt = min(nkxt,n);
    
    for ly = -0.5*(t-1) : 1 : 0.5*(t-1)
        
        % find nearest samples
        nkyt = round(nky+ly);
        
        % compute weighting from kernel
        if( l < 0 ); kwy = max(1-abs(nky-nkyt)*(2/t),0);
        else kwy = kaibes( nky-nkyt, t+1, c, l ); end;
        
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

function ca = kaibes( n, s, c, l )

% --- nearest neighbour interpolation ---
% n = max( min( n, s/2 ), -s/2 - 1/l ); % clip into range [-s/2,s/2)
% ca = c( round( (n+s/2)*l )+1 + 1 ); % nearest neighbour

% --- bilinear interpolation ---
n = max( min( n, s/2 ), -s/2 - 2/l ); % clip into range [-s/2,s/2)
nl = floor(n*l)/l; % push to lower edge of sample interval
inl = floor((n+s/2)*l)+1 + 2; % convert nl to index
ca = (n-nl)*l .* c(inl+1) + (1+(nl-n)*l) .* c(inl); % bilinear

end