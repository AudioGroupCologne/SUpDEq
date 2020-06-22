% [y, not_y, f] = AKfftFilterbank(x, fs, frac, f0)
%
% is an fft-based filterbank with non-overlaping bands of fractional octave
% bandwidth.
%
% Input:
% x         - input signal: mono or multichannel time domain
% fs        - sampling frequency in Hz (default = 44100)
% frac      - fraction of octave determining the bandwith of the filters
%             e.g. frac = 3 -> third octave filters (default = 3)
% f0        - lowest and highest center frequencies to consider
%             (default [32 fs/2])
%
% Output:
% y         - filtered time signal [samples x channels x filterbands]
% not_y     - parts of time signal that are outside the range specified by
%             f0 (first layer: frequencies below)
% f         - 1.-3. col.: ideal center-frequencies, and lower and upper
%                         cut-off frequencies of filterbands
%           - 4.-6. col.: actual frequencies
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
% v.1 05/2015

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
function [y, not_y, f] = AKfftFilterbank(x, fs, frac, f0)

% --------------------------------------------------- check input variables
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('frac', 'var')
    frac = 3;
end
if ~exist('f0', 'var')
    f0 = [32 fs/2];
else
    if length(f0) == 1
        f0(2) = fs/2;
    end
end

if f0(1) == 0
    f0(1) = 32;
    fullrange_low = true;
else
    fullrange_low = false;
end

if f0(2) >= fs/2
    f0(2) = fs/2;
    fullrange_high = true;
else
    fullrange_high = false;
end

N = size(x, 1);
df = fs/N;

% --------------------------------- calculate center and cut-off frequencys

% frequencies according to DIN EN 61260
f = AKfractionalOctaves(frac, f0, fs);
f = f(:,1:3);

% get fft bins of cut-off frequencies
fN = round(f/df)+1;
% make them non-overlapping
fN(:,3) = fN(:,3)-1;
% find bands of zero width
f  = f(fN(:,3)-fN(:,2) >= 0,:,:);
fN = fN(fN(:,3)-fN(:,2) >= 0,:,:);

% calculate actual cut-off frequencies
f(:, 4:6) = (fN-1)*df;


% -------------------------------------- fft -> select wanted parts -> ifft

% get one sided spectrum
X            = fft(x);
[X, is_even] = AKboth2singleSidedSpectrum(X);
N_half       = size(X, 1);

% if f0(2) was fs/2 and spectrum is odd fN(end) is too big
fN(end) = min([fN(end) N_half]);

% get filterbands
Y = zeros(size(X, 1), size(X, 2), length(fN));
for n = 1:length(fN)
    Y(fN(n,2):fN(n,3), :, n) = X(fN(n,2):fN(n,3), :);
end

% save frequencies lower than first filterband
not_Y = zeros(size(X, 1), size(X, 2), 2);
if fN(1) > 1
    not_Y(1:fN(1,2)-1, :, 1) = X(1:fN(1,2)-1, :);
end

% save frequencies higher than last filterband
if fN(end) < N_half
    not_Y(fN(end,3)+1:N_half, :, 2) = X(fN(end,3)+1:N_half, :);
end

% get both sided spectrum
Z     = zeros(N, size(Y,2), size(Y,3));
not_Z = zeros(N, size(Y,2), 2);
for n = 1:size(Y, 2)
    Z(:,n,:)     = AKsingle2bothSidedSpectrum(squeeze(Y(:,n,:)), is_even);
    not_Z(:,n,:) = AKsingle2bothSidedSpectrum(squeeze(not_Y(:,n,:)), is_even);
end

% ifft
y     = ifft(Z,     'symmetric');
not_y = ifft(not_Z, 'symmetric');

if fullrange_low
    y(:,:,1)     = y(:,:,1) + not_y(:,:,1);
    not_y(:,:,1) = 0;
end
if fullrange_high
    y(:,:,end)     = y(:,:,end) + not_y(:,:,end);
    not_y(:,:,end) = 0;
end