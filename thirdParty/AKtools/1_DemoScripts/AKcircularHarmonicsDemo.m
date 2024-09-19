% This Demo script shows some examples of how to use the circular harmonics
% code inside AKtools. The implementation follows [1]
%
% [1] Noam R Shabtai, Gottfried Behler, and Michael Vorlaender: "Generation
%     of a reference radiation pattern of string instruments using
%     automatic excitation and acoustic centering." J. Acoust. Soc. Am.,
%     138(5): EL479-EL456, (2015).
%
% fabian.brinkmann@tu-berlin.de,
% Audio Communication Group, TU Berlin
% 02/2018

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
close all; clear; clc

%% --------------------- 1. compute and plot real valued circular harmonics

% obtain circular harmonics functions at desired angles phi
phi       = 0:359;
[e_n, N] = AKch(2, phi, 'real');

% plot the result
AKf(30,10)
plot(phi, e_n)
xlim([0 phi(end)])
legCell  = cellstr(num2str(N', 'N=%-d'));
legend(legCell, 'location', 'EastOutside')
title  'Real valued circular harmonics for different orders'
xlabel 'azimuth in degree'
ylabel 'amplitude'

%% ------------- 2. Apply circular harmonics transform to a simple function
close all; clear; clc

% set relevant parameters for CH transform and inverse CH transform
phi     = 0:359;
N       = 20;
CHmode  = 'real';

% define a simple function on the circle
h           = ones(size(phi));
h(phi>=180) = -1;

% apply CH transform to obtain CH coefficients d_n
d_n = AKcht(h, false, phi, N, 'complex', 44100, false, CHmode);

% apply inverse CH transform
hCH = AKicht(d_n, false, phi, 'complex', 44100, false, CHmode);
% take real value in case of small imaginary parts caused by numerical noise
hCH = real( hCH );

% plot the result
AKf(20,10)
plot(phi, h)
hold on
plot(phi, hCH)
xlim([0 phi(end)])
legend('before CH transform', 'after CH transform', 'location', 'SouthWest')
title  'Circular function before and after limited order CH transformation'
xlabel 'azimuth in degree'
ylabel 'amplitude'


%% 3. Apply circular harmonics transform to head-related transfer functions
%  This example needs the FABIAN HRTF datase
AKdependencies('Fabian')
close all; clear; clc

% set relevant parameters for CH transform and inverse CH transform
doFFT   = true;
phi     = 0:359;
N       = 35;
CHTmode = 'complex';
CHmode  = 'complex';
compact = false;

% get the HRTFs from the FABIAN database
h = AKhrirInterpolation(phi, 0, 0, 'measured_sh');

% apply CH transform to obtain CH coefficients d_n
d_n = AKcht(h, doFFT, phi, N, CHTmode, 44100, compact, CHmode);

% apply inverse CH transform
hCH = AKicht(d_n, doFFT, phi, CHTmode, 44100, compact, CHmode);
% take real value in case of small imaginary parts caused by numerical noise
hCH = real( hCH );

% plot the result
AKf(20)
subplot(2,1,1)
AKp(h, 'm3d', 'dr', [-20, 20])
subplot(2,1,2)
AKp(hCH, 'm3d', 'dr', [-20, 20])