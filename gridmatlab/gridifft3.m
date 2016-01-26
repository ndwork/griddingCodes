function im = gridifft3( d, k, n, varargin )
% GRIDIIFT3 Compute the three-dimensional inverse-DFT via gridding
%   IM = GRIDIFFT3( D, K, N, [A, L, C, P, T] ) grids the data D to the
%   k-space locations K and returns the s*(NxNxN) matrix that is the
%   inverse-DFT of the gridded (and de-apodised) data, where s is 1 if C is
%   true and A if C is false.  Optional arguments are 
%       A: the oversampling ratio;
%       L: the kernel shape ('tri' - triangle, 'kb' - Kaiser-Bessel);
%       C: a flag indicating whether to crop from the full gridded FOV;
%       P: a flag indicating whether to print reconstruction progress; and
%       T: a flag indicating whether to time each reconstruction operation.
%   By default, a grid over-sampling ratio no smaller than 1.25 and a
%   triangular kernel are used.  Passing the empty matrix [] will indicate
%   the default value for the corresponding argument.
%
%   GRIDIFFT3 uses GRIDMR3D and IFFT3C.
%
% Ethan Johnson, 2013


% Default values
al = 1.25; % grid over-sampling ratio
s = []; % kernel size
a = []; % grid over-sampling ratio specified to gridmr3d
l = 1; % triangle (1: triangle, 2: kaiser-bessel)
thr = 8e-4; % deapodisation threshold
f_prnt = false; % flag: fprintf output of status
f_time = false; % flag: tic-toc timings
f_crop = true; % flag: crop the extended FOV used for gridding


% Parse inputs
if( nargin>3 && ~isempty(varargin{1}) ), al = varargin{1}; end
if( nargin>4 && strcmpi(varargin{2},'kb') ), a = al; l = 2; end
if( nargin>5 && ~isempty(varargin{3}) ), s = varargin{3}; end
if islogical(s), error('Kernel width is now 6th arg.'); end
if( nargin>6 && ~isempty(varargin{4}) ), f_crop = varargin{4}; end
if( nargin>7 && ~isempty(varargin{5}) ), f_prnt = varargin{5}; end
if( nargin>8 && ~isempty(varargin{6}) ), f_time = varargin{6}; end
if numel(n)==1, n = [n n n]; end


if f_prnt, fprintf('Gridding-IFFT reconstruction\n'); end


% Preparation
al = ceil(al*n)/n; aln = round(al*n); % get integer dimensions
a = al(1:length(a)); % update a if appropriate (keep empty if not)
if isempty(s), s = 2*al+1; end % probably too small for good kaiser-bessel
time = ''; % fallback position


% Gridding
if f_time, time=tic; end
im = gridmr3d( d, k, aln, s, a ); clear d;
if f_time, time=toc(time); time=sprintf(': %fs',time); end
if f_prnt, fprintf(' + gridded%s\n',time); end


% IDFT
if f_time, time=tic; end
im = ifft3c( im );
if f_time, time=toc(time); time=sprintf(': %fs',time); end
if f_prnt, fprintf(' + ifft%s\n',time); end


% De-apodise
if f_time, time=tic; end
[x,y,z] = meshgrid( ((-aln(2)/2):(aln(2)/2-1))/aln(2), ...
                    ((-aln(1)/2):(aln(1)/2-1))/aln(1), ...
                    ((-aln(3)/2):(aln(3)/2-1))/aln(3) );
switch l
    case 1
        deap = (sinc((s-1)/2*x).*sinc((s-1)/2*y).*sinc((s-1)/2*z)).^2;
    case 2
        be = pi^2*((s/al*(al-0.5))^2 - 0.8); % (optimal shape parameter)^2
        t = (pi*s*x).^2-be; u = (pi*s*y).^2-be; v = (pi*s*z).^2-be;
        deap = sinc(sqrt(t)/pi) .* sinc(sqrt(u)/pi) .* sinc(sqrt(v)/pi);
end
deap(deap<thr) = thr;
im = im./deap; clear deap z y x;
if f_time, time=toc(time); time=sprintf(': %fs',time); end
if f_prnt, fprintf(' + deapodised%s\n',time); end


% Crop
if f_crop
    if f_time, time=tic; end
    n0 = -floor(n/2)+floor(aln/2);
    cr = { (1:n(2))+n0(2), (1:n(1))+n0(1), (1:n(3))+n0(3) };
    im = im(cr{1},cr{2},cr{3});
    clear cr n0;
    if f_time, time=toc(time); time=sprintf(': %fs',time); end
    if f_prnt, fprintf(' + cropped%s\n',time); end
end


if f_prnt, fprintf('Finished\n'); end


end

