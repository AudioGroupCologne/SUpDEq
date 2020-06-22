% Shows example usages of AKnormalize to normalize data
%
% 11/2016  -  fabian.brinkmann@tu-berlin.de
%
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
close all; clear; clc

% generate some signals
% (we use delayed impulses filtered with a parametric equalizer)
fs = 44100;
x  = AKdirac(1024, 3, [128 256 384]);
x  = AKfilter(x, 'peq', 1000, [20 0 -20], fs, 1, 'hpl', 'tan');

% look at the signals
AKpMulti(x, '2a', 'fs', fs)


%% ------------------- 1. normalization of the time signals to a value of 1

% you can normalize to any value, we use 1 for simplicity
value = 1;

% we normalize the maximum of the time signal, but could also take the
% 'mean', or 'rms'
operation = 'max';

y1 = AKnormalize(x, 'time', operation, 'max',  value);
y2 = AKnormalize(x, 'time', operation, 'min',  value);
y3 = AKnormalize(x, 'time', operation, 'each', value);

AKf
subplot(1,4,1)
    AKp(x, 't2d', 'dr', [-.4 1.6]*value)
    title 'original data'
subplot(1,4,2)
    AKp(y1, 't2d', 'dr', [-.4 1.6]*value)
    title 'normalized largest impulse'
subplot(1,4,3)
    AKp(y2, 't2d', 'dr', [-.4 1.6]*value)
    title 'normalized smallest impulse'
subplot(1,4,4)
    AKp(y3, 't2d', 'dr', [-.4 1.6]*value)
    title 'normalized all impulses'


%% ----------- 2. normalization of the broad band magnitude spectra to 0 dB

% you can normalize to any value, we use 1 for simplicity
value = 0;

% we normalize the maximum of the time signal, but could also take the
% 'mean', or 'rms'
operation = 'max';

y1 = AKnormalize(x, 'dB', operation, 'max',  value);
y2 = AKnormalize(x, 'dB', operation, 'min',  value);
y3 = AKnormalize(x, 'dB', operation, 'each', value);

AKf
subplot(1,4,1)
    AKp(x, 'm2d', 'dr', [-45 25])
    title 'original data'
subplot(1,4,2)
    AKp(y1, 'm2d', 'dr', [-45 25])
    title 'normalized largest magnitude'
subplot(1,4,3)
    AKp(y2, 'm2d', 'dr', [-45 25])
    title 'normalized smallest magnitude'
subplot(1,4,4)
    AKp(y3, 'm2d', 'dr', [-45 25])
    title 'normalized all magnitudes'


%% ---------- 3. normalization of the narrow band magnitude spectra to 0 dB

% you can specify a frequency range for normalization either by center
% frequency and width in fractional octaves, or by upper and lower
% frequency bounds. We use the first here
f    = 1000;  % center frequ., you could also try [1000 10000]
frac = 1/3;   % width is ignored if f = [a b]

% you can normalize to any value, we use 1 for simplicity
value = 0;

% we normalize the maximum of the time signal, but could also take the
% 'mean', or 'rms'
operation = 'mean';

AKnormalize(x, 'dB', operation, 'each', value, f, frac, fs, true);
