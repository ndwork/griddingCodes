function gb = GAMMAB(varargin)
% GAMMAB Return the value of \bar{gamma}=gamma/2/pi.
%   GAMMAB returns the value in Hz/T
%   
%   GAMMAB(UNT) returns the value in units of UNT (case insensitive).
%       supported: 'Hz/T' [default], 'Hz/G', 'rad/T', 'rad/G'
%   
%   GAMMAB(UNTF,UNTM) returns the value in units of UNTF/UNTM.
%       supported UNTF: 'Hz' [default], 'kHz', 'MHz', 'rad', 'krad', 'Mrad'
%       supported UNTM: 'T', 'G', 'mT', 'uT'
%
%   N.b. 'rad' refers to radians (per second).

s = 1;
if( nargin > 1 )
    switch lower(varargin{1}),
        case 'khz', s = s * 1e-3;% display('KILOHERTZ!');
        case 'mhz', s = s * 1e-6;% display('MEGAHERTZ!');
        case 'rad', s = s * 2*pi;% display('RADIANS!');
        case 'krad', s = s * 2*pi*1e-3;% display('KILORADIANS!');
        case 'mrad', s = s * 2*pi*1e-6;% display('MEGARADIANS!');
        case 'hz', s = s * 1;
        otherwise, fprintf('***Warning, unrecognised units for GAMMAB\n');
    end
    switch lower(varargin{2}),
        case 'g', s = s * 1e-4;% display('GAUSS!');
        case 'mt', s = s * 1e-3;% display('MILLITESLA!');
        case 'ut', s = s * 1e-6;% display('MICROTESLA!');
        case 't', s = s * 1;
    end
elseif( nargin > 0 )
    switch lower(varargin{1})
        case 'hz/g', s = s * 1e-4;
        case 'rad/t', s = s * 2*pi;
        case 'rad/g', s = s * 2*pi*1e-4;
        case 'hz/t', s = s * 1;
        otherwise, fprintf('***Warning, unrecognised units for GAMMAB\n');
    end
end

gb = s * 42.576e6; % [Hz/T]

end