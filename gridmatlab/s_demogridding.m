%% s_demogridding
%  Ethan Johnson, 2015

%% Gradient shape

% --- Set up ---

% Simulation parameters
dt = 1e-3; % time step [ms]
pf = 64; % grid padding size
t2 = 0.13; % T2 [ms]
% -- collection 1 --------------------------------------------------------
al = 1.25; % gridding oversampling factor
sh = 'tri'; % gridding kernel shape
s = []; % kernel width (ignored if triangle)
% % ------------------------------------------------------------------------
% % -- collection 4 --------------------------------------------------------
% al = 2; % gridding oversampling factor
% sh = 'tri'; % gridding kernel shape
% s = []; % kernel width (ignored if triangle)
% ------------------------------------------------------------------------
% % -- collection 2 --------------------------------------------------------
% al = 1.25; % gridding oversampling factor
% sh = 'kb'; % gridding kernel shape
% s = 6; % kernel width (ignored if triangle)
% % ------------------------------------------------------------------------
% % -- collection 3 --------------------------------------------------------
% al = 1.5; % gridding oversampling factor
% sh = 'kb'; % gridding kernel shape
% s = 5; % kernel width (ignored if triangle)
% % ------------------------------------------------------------------------

% Phys. / acq. parameters
fov = 32e-2; % field of view [m]
rf = 0.15/2; % field of view undersampling factor
sm = 200; % slew-max [mT/m/ms]
gm = 33; % g-max [mT/m]
tro = 1; % read-out time [ms]
rscc = 1.59043899; % resolution scaling (FWHM) for corner-cut k-space (3D)

% Display options
% --



% --- Make gradient ---
% This gradient is for a radial trajectory

% Variable-density / SNR-optimised gradients (Nagel '09)
trmp = ceil( gm/sm /dt)*dt; % ramp-up time [ms] (guess: slewing max rate)
tflt = round( 130e-3 /dt)*dt; % platea time [ms] (inferred Jo12fig1)
t0 = trmp+tflt; k0 = GAMMAB('kHz','mT')*gm*(0.5*trmp + tflt);
t = t0+dt:dt:tro;
g = [ linspace(0,gm,trmp/dt) ...
      gm*ones(1,tflt/dt) ...
      k0^2 * gm * ( 3*GAMMAB('kHz','mT')*k0^2*gm*(t-t0) + k0^3 ).^(-2/3) ];
clear trmp tflt t0 t;
ag = sum(g*dt); % gradient area [ms mT/m]



% --- Ascertain prameters ---

% Figure phys. params.
dx = 1 / ( GAMMAB('kHz','mT')*ag * 1e-3 ); % resolution [mm]
bwm = GAMMAB('kHz','mT') * max(g(:)) * fov; % max bandwidth [kHz]
nx = round((fov*1e3)/dx); nnn = nx*[1 1 1];

% Figure acq. params.
nspo = ceil(pi*(nx*rf)^2); % for 3d: 4*pi*r^2 = 4*pi*(n/2)^2 = pi*n^2

% Other shortcuts
[ng,nsa] = size(g); % nbr of gradient waveforms, nbr of samples per wvfrm
n2 = numel(t2); % nbr of t2 values (add sth. here to condition shape of t2s ?)
n0 = floor(nx/2)+1; n0p = floor(nx*pf/2)+1; % indices for central sample



% --- Make k-space trajectories ---

% Make kr
kr = cumsum(g,2); kr = permute(bsxfun(@times,kr,1./kr(:,end)),[2 3 4 1]);

% Make ke
n = (0:nspo-1).';
phs = -acos(2*n/nspo - 1 + 1/nspo);
ths = mod(1.285*sqrt(2*nspo)*asin(2*n/nspo - 1 + 1/nspo),2*pi);
rhs = 0.5*ones(size(phs)); dc = ones(size(rhs));
clear n;
ke = [rhs.*sin(phs).*cos(ths) rhs.*sin(phs).*sin(ths) rhs.*cos(phs)];

% Assemble k-space points from the end-points and radial trajectory
k = bsxfun(@times, kr, permute(ke,[3 1 2]));

% Compute weights
wr = [1/8*kr(2,:,:,:).^3 ; diff( kr.^3 )]; % density (spherical shells)
w = bsxfun(@times,dc.',wr);



% --- Impulse response ---

% % Make up data
t = (0:dt:tro-dt).'; e2 = exp(-t/t2);
d = e2 * ones(1,nspo);
clear t e2;

thm = 0.5; tom = 0.01; % full-width thresholds

% Reconstruct
fprintf('Reconstructing impulse response\n---->\n'); time=tic;
    
% Reconstruct impulse response
timei = tic;
if strcmpi(sh,'tri'), s=2*al+1; end
im = gridifft3(d.*w,k,nnn,al,sh,s,true,false,false);
  % d is data
  % w is density compensation
  % k is the k space trajectory
  % nnn is the image size
  % al is the oversampling ratio
  % sh is the shape of the kernel
  % s is the width of the kernel
timei = toc(timei);
    
% Extract and sinc-interpolate cross-sections
timej = tic;
cxp = ifftc(padarray(fftc(im(n0,:,n0)),[0 (pf-1)*nnn(2)/2 0],0,'both'));
cyp = ifftc(padarray(fftc(im(:,n0,n0)),[(pf-1)*nnn(1)/2 0 0],0,'both'));
czp = ifftc(padarray(fftc(im(n0,n0,:)),[0 0 (pf-1)*nnn(1)/2],0,'both'));
timej = toc(timej);
    
% Ascertain FW_M
timek = tic;
hx = find(abs(cxp(n0p:end))<=thm*cxp(n0p),1); if isempty(hx), hx=nan; end
hy = find(abs(cyp(n0p:end))<=thm*cyp(n0p),1); if isempty(hy), hy=nan; end
hz = find(abs(czp(n0p:end))<=thm*czp(n0p),1); if isempty(hz), hz=nan; end
fwhm = 2*[hx;hy;hz]; % figure in terms of half-widths and double
hx = find(abs(cxp(n0p:end))<=tom*cxp(n0p),1); if isempty(hx), hx=nan; end
hy = find(abs(cyp(n0p:end))<=tom*cyp(n0p),1); if isempty(hy), hy=nan; end
hz = find(abs(czp(n0p:end))<=tom*czp(n0p),1); if isempty(hz), hz=nan; end
fwom = 2*[hx;hy;hz]; % figure in terms of half-widths and double
timek = toc(timek);

fprintf('   + %f (grid)|%fs (interp.)|%fs (meas.)\n',timei,timej,timek);

fprintf('<----\n');
time=toc(time); fprintf(' + %fs\n',time); clear time;

% Rescale FW_M to physical units
fwhm = fwhm * (dx/pf); fwom = fwom * (dx/pf);



% --- Prepare for display ---

% Rescale cross sections and images to have unit peak magnitude
im = bsxfun(@times,im,1./im(n0,n0,n0));
cxp = bsxfun(@times,cxp(:),1./cxp(n0p));
% -- not going to display cy or cz
% Extract un-interpolated cross-sections
cx = reshape(im(n0,:,n0),[],1);
% Make spatial grid
x = linspace(-fov/2,fov/2-dx*1e-3,nx) *1e3; % spatial scale
xp = linspace(-fov/2,fov/2-dx/pf*1e-3,nx*pf) *1e3; % interp'd spatial scale



% --- Display ---

npx = 7; % width of FW_M demarcation

% Cross-sections
fig=figure;
yl = 10.^[-4 0];
% (--- full FOV ---)
subplot(211);
semilogy(x(n0+npx*[-1 1]),thm*[1 1],'k','LineWidth',0.5); % half-max line
hold on;
semilogy(x(n0+npx*[-1 1]),tom*[1 1],'k','LineWidth',0.5); % 1%-max line
semilogy(x(n0-npx)*[1 1],yl,'k:', x(n0+npx)*[1 1],yl, 'k:','LineWidth',0.5); % lines delimiters
for i=1:ng, semilogy( x,abs(cx(:,i)),'.' ); semilogy( xp,abs(cxp(:,i)),'LineWidth',1 ); end; clear i;
hold off;
xlabel('{\itx} [mm]'); ylabel('ampl.'); xlim(x([1 end])); ylim(yl);
% (--- zoom FOV ---)
subplot(212);
nzm = n0 + (-npx:npx); nzmp = n0p + (-npx*pf:npx*pf);
semilogy(x(n0+npx)*[0 1],thm*[1 1],'k','LineWidth',0.5); % half-max line
hold on;
semilogy(x(n0+npx)*[0 1],tom*[1 1],'k','LineWidth',0.5); % 1%-max line
% - Reference impulse responses (sinc / 'sphinc') -
clrsinc = 0.92* 0.6/0.92* [255 255 255]/256; clrsphinc = 0.7* 0.4/0.7* [255 255 255]/256;
px = xp(nzmp)/dx; % spatial points in units of pixels
plot( xp(nzmp), abs( sinc(px) ), '-','Color',clrsinc); % sinc (cartesian impulse response)
[~,ixz] = min(abs(xp-(x(n0)-dx))); pxz = xp(ixz)/dx;
text(xp(ixz),max(abs(sinc(pxz)),min(yl)),'sinc','Color',clrsinc); % label sinc
clear ixz pxz;
plot( xp(nzmp), abs( 3./(pi*px).^2 .* (sinc(px)-cos(pi*px)) ), '-','Color',clrsphinc); % sphinc (3d corner-cut impulse response)
[~,ixz] = min(abs(xp-(x(n0)-dx*1.43029668))); pxz = xp(ixz)/dx;
text(xp(ixz),max( abs(3./(pi*pxz).^2 .* (sinc(pxz)-cos(pi*pxz))) ,min(yl)),'{\ith}','Color',clrsphinc); % label sphinc
clear ixz pxz clrsinc clrsphinc px;
% -
% - Impulse responses -
for i=1:ng, semilogy( x(nzm),abs(cx(nzm,i)),'.' ); semilogy( xp(nzmp),abs(cxp(nzmp,i)) ); end; clear i;
% - 
hold off;
xlabel('{\itx} [mm]'); ylabel('ampl.'); xlim(x(n0+npx*[-1 1])); ylim(yl);
clear nzm nzmp;
clear yl;
annotation('textbox',[0 0 1 1],... % title
           'String',sprintf('gridding: {\\it\\alpha}=%.2f, {\\its}=%1f, shape=%s',al,s,sh),...
           'HorizontalAlign','center','VerticalAlign','bottom',...
           'EdgeColor','none',...
           'FontName',get(0,'DefaultAxesFontName'),'FontSize',get(0,'DefaultAxesFontSize'));
pause; try close(fig); catch; end; clear fig;

clear npx;
