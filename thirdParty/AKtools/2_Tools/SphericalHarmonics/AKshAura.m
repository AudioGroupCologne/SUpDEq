% [y, azI, elI] = AKshAura(x, az, el, fnm1, fnm2, SHTmode, doIFFT, isEven, N, toa1, toa2)
% renders an auralization for the trajectory defined by az and el based on
% the input audio file x that is convolved with impulse responses obtained
% from spherical harmonics (SH) coefficients fnm1 and fnm2.
%
% See AKsphericalHarmonicsTransformDemo.m for examples
%
% auralization of horizontal plane
% AKshAura(rand(5*44100,1)*2-1, [0 360], [0 0], f_nm1, f_nm2);
%
% auralization of median plane
% AKshAura(rand(5*44100,1)*2-1, [0 0 NaN 180 180], [-90 90 NaN 90 -90], f_nm1, f_nm2);
%
% auralization of frontal plane
% AKshAura(rand(5*44100,1)*2-1, [90 90 NaN 270 270], [-90 90 NaN 90 -90], f_nm1, f_nm2);
%
% INPUT:
% x       - One, or two channel time domain audio signal.
%           One channel= one column
% az      - Azimuth in degree for auralization according to SOFA coordinate
%           convention (See AKsh.m). Scalar or vector.
%           NaNs define segments in az (see example calls above).
% el      - Elevation in degree for auralization according to SOFA
%           coordinate convention. Scalar or vector of size(az).
% fnm1    - fnm spherical harmonics coefficients of first channel IRs.
%           See AKsh.m for help
% fnm2    - fnm spherical harmonics coefficients of second channel IRs, or
%           empty array (default = []).
% SHTmode - See AKsh.m for help (default = 'db_unwrap').
% doIFFT  - See AKish.m for help (default = true).
% isEven  - See AKsh.m for help (default = 1).
% N       - Block size for convolution (default = 512).
% toa1    - Spherical harmonics coefficients describing the left ear time
%           of arrival in samples. If this is passed, the time of arrival
%           is obtained from the coefficients, and the IRs obtained from
%           fnm1 are delayed accordingly
%           (Format as in HUTUBS HRTF database)
% toa2    - Spherical harmonics coefficients describing the right ear time
%           of arrival.
%
% OUTPUT:
% y       - Input signal x convolved with IRs from fnm coefficients
% azI     - actual aziumth angles used for convolution
% elI     - actual elevation angles used for convolution
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
% 5/2016

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License.
function [y, azI, elI] = AKshAura(...
    x, az, el, fnm1, fnm2, SHTmode, doIFFT, isEven, N, toa1, toa2)

if ~exist('N', 'var')
    N = 512;
end
if ~exist('isEven', 'var')
    isEven = true;
end
if ~exist('doIFFT', 'var')
    doIFFT = true;
end
if ~exist('SHTmode', 'var')
    SHTmode = 'db_unwrap';
end
if ~exist('fnm2', 'var')
    fnm2 = [];
end
if ~exist('toa1', 'var')
    toa1 = [];
end
if ~exist('toa2', 'var')
    toa2 = [];
end

if numel(az) ~= numel(el)
    error('az and el must have the same size')
end

if numel(az) > 1
    
    % get right format for az and el
    az = reshape(az, [numel(az), 1]);
    el = reshape(el, [numel(el), 1]);
    
    % find segments in azimuth and elevation
    Sid = find(isnan(az)==1); numel(az);
    S   = numel(Sid) + 1;
    if any(~isnan(el(Sid)))
        error('NaN must be at the same positions in az and el.')
    end
    Sid = [0; Sid; numel(az)+1];
    
    % get number of Blocks for convolution
    L = ceil(size(x,1)/N);
    if mod(L, S)
        L = L + S-mod(L, S);
    end
    
    % zero Pad input signal to integer divisor of block size N
    x(end+1:N*L,:) = 0;
    
    % Interpolate trajectory
    azI = zeros(L,1);
    elI = azI;
    for ss = 1:S
        azI((ss-1)*L/S+1:ss*L/S) = interp1(az(Sid(ss)+1:Sid(ss+1)-1), linspace(1, numel(Sid(ss)+1:Sid(ss+1)-1), L/S));
        elI((ss-1)*L/S+1:ss*L/S) = interp1(el(Sid(ss)+1:Sid(ss+1)-1), linspace(1, numel(Sid(ss)+1:Sid(ss+1)-1), L/S));
    end
    
    % get IRs for first channel
    h1 = get_IRs(fnm1, toa1, azI, elI, doIFFT, SHTmode, isEven);
    
    % length of IRs
    Nh = size(h1,1);
    
    % get irs of second channel
    if ~isempty(fnm2)
        h2 = get_IRs(fnm2, toa2, azI, elI, doIFFT, SHTmode, isEven);
        
        % match size of input signal
        if size(x,2) == 1
            x = repmat(x, [1 2]);
        else
            x = x(1:end, 1:2);
        end
    else
        x = x(:,1);
    end
    
    
    % allocate output signal
    y = zeros(N*L+Nh-1, size(x,2));
    
    % fft filt
    if size(x,2) == 1
        for ll = 1:L
            y((ll-1)*N+1:ll*N+Nh-1)   = y((ll-1)*N+1:ll*N+Nh-1)   + fftfilt(x((ll-1)*N+1:ll*N),    [h1(:,ll); zeros(N-1, 1)]);
        end
    else
        for ll = 1:L
            y((ll-1)*N+1:ll*N+Nh-1,:) = y((ll-1)*N+1:ll*N+Nh-1,:) + fftfilt(x((ll-1)*N+1:ll*N,:), [[h1(:,ll) h2(:,ll)]; zeros(N-1, 2)]);
        end
    end   
    
else
    
    % get IRs for first channel
    h1 = get_IRs(fnm1, toa1, az, el, doIFFT, SHTmode, isEven);
    
    % get irs of second channel
    if ~isempty(fnm2)
        h2 = get_IRs(fnm2, toa2, az, el, doIFFT, SHTmode, isEven);
        
        % match size of input signal
        if size(x,2) == 1
            x = repmat(x, [1 2]);
        else
            x = x(1:end, 1:2);
        end
    else
        x = x(:,1);
    end
    
    % fft filt
    if size(x,2) == 1
        y = fftfilt(h1,      [x; zeros(size(h1,1),1)]);
    else
        y = fftfilt([h1 h2], [x; zeros(size(h1,1),2)]);
    end   
    
    azI = az;
    elI = el;

end

% normalize to match input level
y = y * max(abs(x(:))) / max(abs(y(:)));

end

function h = get_IRs(fnm, toa, az, el, doIFFT, SHTmode, isEven)
% get IRs for first channel
h = AKisht(fnm, doIFFT, [az 90-el], SHTmode, isEven);
    
if ~isempty(toa)
    d = AKisht(toa, false, [az 90-el], ...
        'complex', true, false, 'real');
    h = AKfractionalDelayCyclic(h, d);
end
end
