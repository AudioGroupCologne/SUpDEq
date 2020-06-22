% [h, h_inv] = AKsphericalHeadDiffuse(a, phase, N, nSamples, fs, c)
%
% calculates the diffuse field impulse response, and incerted diffuse field
% impulse response of a spherical head model according to eq. (11) in [1].
%
% e.g.
% [h, h_inv] = AKsphericalHeadDiffuse;
% AKp([h h_inv])
% calucalates and plots the impulse responses using the default parameters
%
%
% I N P U T
% a        - radius of the spherical head in meter (default = 0.0875)
% phase    - desired phase behaviour
%            'min' - generate minimum phase filter (default)
%            'lin' - generate linear phase filter
% N        - Spherical Harmonics order (default = 50)
% nSamples - length of impulse response in samples (default = 128)
% fs       - sampling rate in Hz (default = 44100)
% c        - speed of sound in m/s (default = 343)
%
%
% O U T P U T
% h        - time domain diffuse field response of the spherical head
% h_inv    - inverted version of h
%
%
% [1] Zamir Ben-Hur, Fabian Brinkmann, Jonathan Sheaffer, Stefan Weinzierl
%     and Boaz Rafaely: "Spectral equalization in binaural signals
%     represented by order-truncated spherical harmonics." J. Acoust.Soc.
%     Am., 141(6):4087-4096, 2017.
%
%
% 01/2018 - fabian.brinkmann@tu-berlin.de

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
function [h, h_inv] = AKsphericalHeadDiffuse(a, phase, N, nSamples, fs, c)


if ~exist('N', 'var')
    N = 50;
end
if ~exist('a', 'var')
    a = .0875;
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

% ----------------------------------------------- generate complex spectrum
% frequencies and wave numbers to be calculated (do not calculate 0 Hz)
f  = (fs/nSamples:fs/nSamples:fs/2)';
k  = 2*pi*f/c;
kr = k*a;

% eq. (9) from [1]
nn   = 0:N;

j_n  = AKshRadial(kr, 'bessel', 1, nn, false); % spherical Bessel function of first kind
j_nd = AKshRadial(kr, 'bessel', 1, nn, true);  % derived spherical Bessel function of first kind
h_n  = AKshRadial(kr, 'hankel', 1, nn, false); % spherical Hankel function of first kind
h_nd = AKshRadial(kr, 'hankel', 1, nn, true);  % derived spherical Hankel function of first kind

b_n = repmat(4*pi*1i.^nn, [numel(kr) 1]) .* ( j_n - j_nd./h_nd .* h_n );

% eq. (11)
b_n_sq  = abs(b_n).^2;

nn   = 2*(0:N)+1;
nn_b = repmat(nn, [numel(kr) 1]) .* b_n_sq(:,1:N+1);
H    = 1/(4*pi) * sqrt( sum( nn_b, 2 ) );

% ----------------------------------------------- generate impulse response
% add 0 Hz bin
H = [1; H];

% make sure bin at fs/2 is real
if f(end) == fs/2
    H(end) = abs(H(end));
end

% mirror the spectrum
H     = AKsingle2bothSidedSpectrum( H,    1-mod(nSamples, 2) );

% get zero phase impulse response
h     = ifft(H,    'symmetric');
h_inv = ifft(1./H, 'symmetric');

% generate desired phase
if strcmpi(phase, 'lin')
    h     = AKphaseManipulation(h,     fs, phase, 0, false);
    h_inv = AKphaseManipulation(h_inv, fs, phase, 0, false);
else
    err = 2;
    NFFTdouble = 0;
    while db(err(1)) > 0.01
        NFFTdouble = NFFTdouble + 1;
        [h, err] = AKphaseManipulation(h,     fs, phase, NFFTdouble, false);
         h_inv   = AKphaseManipulation(h_inv, fs, phase, NFFTdouble, false);
    end
end