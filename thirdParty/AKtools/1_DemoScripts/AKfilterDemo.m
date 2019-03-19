% Demonstration for filtering audio signals with AKfilter.m using filter
% functions from Frank Schultz
%
% A python implementation of most filters can be found here:
% http://nbviewer.jupyter.org/github/sonible/nb/blob/master/iir_filter/index.ipynb
%
% 3/2015 initial dev. - fabian.brinkmann@tu-berlin.de
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

% suppress unused variables warning for demo script
%#ok<*ASGLU>

% Generate input signal
x = AKdirac(2^11);

%% low pass / high pass
%  all input arguments exept 'type' can be scalars or vectors
AKf

y = AKfilter(x, 'lp', 1000, 0, 44100, 4, {'LR' 'butter'});
subplot(2,1,1)
AKp(y, 'm2d', 'c', 'cyc')
title('Linkwitz-Riley and butterworth low passes')

y = AKfilter(x, 'hp', [125 250 500 1000 2000 4000 8000], 0, 44100, 4, 'butter');
subplot(2,1,2)
AKp(y, 'm2d', 'c', 'cyc')
title('High passes at different cut-off frequencies')

%% band pass / band stop
AKf

y = AKfilter(x, 'bp', [500 2000], 0, 44100, 16, 'butter');
subplot(2,1,1)
AKp(y, 'm2d', 'c', 'cyc')
title('16th order butterworth band pass')

y = AKfilter(x, 'bs', [500 2000], 0, 44100, 16, 'butter');
subplot(2,1,2)
AKp(y, 'm2d', 'c', 'cyc')
title('16th order butterworth band stop')

%% cross over

y = AKfilter(x, 'xover', 1000,0, 44100, 8, 'LR');
AKp(y(:,:,1), 'm2d', 'c', 'b')
AKp(y(:,:,2), 'm2d', 'c', 'r')
AKp(y(:,:,1)+y(:,:,2), 'm2d', 'ls', '--')
title('8th order Linkwitz Riley cross over network')
legend('lowpass', 'highpass', 'lp + hp', 'location', 'SouthEast')

%% Parametric equalizers
AKf

y = AKfilter(x, 'peq', 1000, -9:3:9, 44100, 4, 'hpl', 'tan');

subplot(3,1,1)
AKp(y, 'm2d', 'c', 'cyc', 'dr', [-10 10], 'x', [100 10000])
title('Half pad loss PEQ')

y = AKfilter(x, 'peq', 1000, -9:3:9, 44100, 4, 'sym', 'tan');

subplot(3,1,2)
AKp(y, 'm2d', 'c', 'cyc', 'dr', [-10 10], 'x', [100 10000])
title('Symmetrical PEQ')

y = AKfilter(x, 'peq', 1000, -9:3:9, 44100, 4, 'con', 'tan');

subplot(3,1,3)
AKp(y, 'm2d', 'c', 'cyc', 'dr', [-10 10], 'x', [100 10000])
title('Constant Q PEQ')

%% LowShelve / highshelve
AKf

y = AKfilter(x, 'ls', 1000, -9:3:9, 44100, 2, 'low');

subplot(3,1,1)
AKp(y, 'm2d', 'c', 'cyc', 'dr', [-10 10])
title('2nd order lowshelves (type ''low'')')

y = AKfilter(x, 'hs', 1000, -9:3:9, 44100, 1, 'low');

subplot(3,1,2)
AKp(y, 'm2d', 'c', 'cyc', 'dr', [-10 10])
title('2nd order highshelves (type ''low'')')

y = AKfilter(x, 'ls', 1000, [9 9 9 -9 -9 -9], 44100, 2, {'low' 'mid' 'high' 'low' 'mid' 'high'});

subplot(3,1,3)
AKp(y, 'm2d', 'c', 'bkrbkr', 'dr', [-10 10])
title('2nd order lowshelves (all types)')
legend('low', 'mid', 'high', 'location', 'SouthEast')

%% Fractional octave filter bank (far from perfectly reconstructing!)
AKf

[y, f] = AKfilter([x x], 'fracOct', [250 8000], 0, 44100, [4 8], 'Butter', 3);

subplot(2,1,1)
AKp(squeeze(y(:,1,:)), 'm2d')
title('4th order Butterworth 3rd Octave filterbank')

subplot(2,1,2)
AKp(squeeze(y(:,2,:)), 'm2d')
title('8th order Butterworth 3rd Octave filterbank')

%% FFT brickwall fractional octave filter bank

[y, f] = AKfilter(AKdirac(2^13), 'FFTfracOct', [0 22050], 0, 44100, 3);

AKf
subplot(2,1,1)
AKp(squeeze(y), 'm2d', 'c', 'cyc')
title('3rd octave FFT brickwall filter bank')

subplot(2,1,2)
AKp(sum(squeeze(y), 2), 'm2d')
title('3rd octave FFT brickwall filter bank (sum)')

%% Peak/Notch filters from experiments by Moore1989
%  (see AKfilter for exact reference)

y = AKfilter(x, 'moore1989', 1000, -9:3:9, 44100, .5);
AKp(y, 'm2d', 'c', 'cyc', 'dr', [-10 10])

y = AKfilter(x, 'moore1989', 8000, -9:3:9, 44100, .5);
AKp(y, 'm2d', 'c', 'cyc', 'dr', [-10 10], 'x', [100 20000])

title('') % title
