% [g_SHF, w_taper, ww_taper] = AKshEqualization(N, N_high, r0, phase, nSamples, fs, c, taper)
%
% calculates an FIR equalization filter to acount for order truncation in
% spherical harmonics based binaural signals according to eq. (12) in [1].
%
% e.g.
% g_SHF = AKshEqualization(12, 30);
% AKp(g_SHF)
% calucalates and plots a filter for equalizing order 12 to 30 (cf. Fig. 3
% in [1])
%
% I N P U T
% N        - low truncation order of binaural signals (Integer >= 0)
% N_high   - high truncation order of binaural signals (Integer > N) 
% r0       - radius of the spherical head in meter (default = 0.0875)
% phase    - desired phase behaviour of the equalization filter g_SHF.
%            'min' - generate minimum phase filter (default)
%            'lin' - generate linear phase filter
% nSamples - length of impulse response in samples (default = 128)
% fs       - sampling rate in Hz (default = 44100)
% c        - speed of sound in m/s (default = 343)
% taper    - apply tapering according to [2]. (default = false)
%
% O U T P U T
% g_SHF    - time domain equalization filter
% w_taper  - tapering window of size [N+1 x 1] according to [2]. If
%            taper=false, w_taper only holds ones.
% ww_taper - tapering window of size [(N+1)^2 x 1] according to [2]. This
%            is the same as w_taper, but the entries for each order are
%            repeated, i.e., the w_taper(2), which holds the entry for the
%            first order is repeated three times for degrees -1, 0, and 1,
%            w_taper(3) is repeated five times, etc.
%
% [1] Zamir Ben-Hur, Fabian Brinkmann, Jonathan Sheaffer, Stefan Weinzierl
%     and Boaz Rafaely: "Spectral equalization in binaural signals
%     represented by order-truncated spherical harmonics." J. Acoust.Soc.
%     Am., 141(6):4087-4096, 2017.
% [2] Christoph Hold, Hannes Gamper, Ville Pulkki, Nikunj Raghuvanshi, and
%     Ivan J. Tashev (2019) "Improving Binaural Ambisonics Decoding by
%     Spherical Harmonics Domain Tapering and Coloration Compensation." In
%     IEEE Int. Conf. Acoustics, Speech and Signal Processing (ICASSP),
%     Brighton, UK, 2019 (pp. 261-265). 

%
% 10/2017 - v1 fabian.brinkmann@tu-berlin.de
% 04/2020 - v2 fabian.brinkmann@tu-berlin.de - added tapering

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [g_SHF, w_taper, ww_taper] = AKshEqualization(N, N_high, r0, phase, nSamples, fs, c, taper)

if ~exist('r0', 'var')
    r0 = .0875;
end
if ~exist('phase', 'var')
    phase = 'min';
end
if ~exist('nSamples', 'var')
    nSamples = 128;
end
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('c', 'var')
    c = 343;
end
if ~exist('taper', 'var')
    taper = false;
end

% ----------------------------------------------------- get tapering window
w_taper  = ones((N+1), 1);
ww_taper = ones((N+1)^2, 1);

if taper
    % window with one entry per order
    w_hann  = hann(2 * (floor( (N+1)/2) ) + 1);
    w_taper(end-(floor( (N-1)/2) ):end) = w_hann(end-(floor( (N+1)/2) ):end-1);
    
    % window with one entry per order and degree
    for nn = 0:N
        MM = max(1, nn^2+1) : (nn+1)^2;
        ww_taper(MM) = w_taper(nn+1);
    end
end

% ----------------------------------------------- generate complex spectrum
% frequencies and wave numbers to be calculated (do not calculate 0 Hz)
f  = (fs/nSamples:fs/nSamples:fs/2)';
k  = 2*pi*f/c;
kr = k*r0;

% eq. (9) from [1]
nn   = 0:N_high;

j_n  = AKshRadial(kr, 'bessel', 1, nn, false); % spherical Bessel function of first kind
j_nd = AKshRadial(kr, 'bessel', 1, nn, true);  % derived spherical Bessel function of first kind
h_n  = AKshRadial(kr, 'hankel', 1, nn, false); % spherical Hankel function of first kind
h_nd = AKshRadial(kr, 'hankel', 1, nn, true);  % derived spherical Hankel function of first kind

b_n = repmat(4*pi*1i.^nn, [numel(kr) 1]) .* ( j_n - j_nd./h_nd .* h_n );

% eq. (11)
b_n_sq  = abs(b_n).^2;

nn   = 2*(0:N)+1;
nn_b = repmat(nn, [numel(kr) 1]) .* b_n_sq(:,1:N+1);
p_N  = 1/(4*pi) * sqrt( sum( nn_b .* w_taper', 2 ) );

nn       = 2*(0:N_high)+1;
nn_b     = repmat(nn, [numel(kr) 1]) .* b_n_sq(:,1:N_high+1);
p_N_high = 1/(4*pi) * sqrt( sum( nn_b, 2 ) );

% eq. (12)
G_SHF = p_N_high ./ p_N;

% ----------------------------------------------- generate impulse response
% add 0 Hz bin
G_SHF = [1; G_SHF];

% mirror the spectrum
G_SHF = AKsingle2bothSidedSpectrum( G_SHF, 1-mod(nSamples, 2) );

% get zero phase impulse response
g_SHF = ifft(G_SHF, 'symmetric');

% generate desired phase
if strcmpi(phase, 'lin')
    g_SHF = AKphaseManipulation(g_SHF, fs, phase, 0, false);
else
    err = 2;
    NFFTdouble = 0;
    while db(err(1)) > 0.01
        NFFTdouble = NFFTdouble + 1;
        [g_SHF, err] = AKphaseManipulation(g_SHF, fs, phase, NFFTdouble, false);
    end
end