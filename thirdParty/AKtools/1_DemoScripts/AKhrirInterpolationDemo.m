% Example for HRIR interpolation on the AK FABIAN HRIR data [1], [2] set
% for multiple head-above-torso orientations (HATO)
%
% The source position is interpolated using spherical harmonics
% coefficients. The HATO is interpolated using inverse distance weighting,
% as described in [3] using the magnitude and unwrapped phase spectra.
%
% [1] JASA???
% [2] Deposit Once
% [3] F Brinkmann, R Roden,A Lindau, S Weinzierl "Audibility and inter-
%     polation of head-above-torso orientation in binaural technology."
%     IEEE J. Sel. Topics Signal Process., 9(5), 931-942, (2015). 

%
% 6/2016 - fabian.brinkmann@tu-berlin.de

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

%% ---------------------------------------- check if data base is available
AKdependencies('FABIAN')


%% ------------------------------------------------------------ single HRIR

% interpolation is done inside this function.
[l, r] = AKhrirInterpolation(90, 0, 0);

% quick plot
AKf
subplot(2,1,1)
    AKp([l r], 't2d', 'c', [.2 .2 .8; .8 .2 .2])
    legend('left ear', 'right ear', 'location', 'NorthEast')
subplot(2,1,2)
    AKp([l r], 'm2d', 'c', [.2 .2 .8; .8 .2 .2], 'N', 44100/20)

% listen via headphones if you want
% soundsc(fftfilt([l r], randn(44100,1)), 44100)


%% ------------------------------------- single directional impuse response

% interpolation is done inside this function.
[l, r] = AKhrirInterpolation(90, 0, 0, 'dir');

% quick plot
AKf
subplot(2,1,1)
    AKp([l r], 't2d', 'c', [.2 .2 .8; .8 .2 .2])
    legend('left ear', 'right ear', 'location', 'NorthEast')
subplot(2,1,2)
    AKp([l r], 'm2d', 'c', [.2 .2 .8; .8 .2 .2], 'N', 44100/20)

% listen via headphones if you want
% soundsc(fftfilt([l r], randn(44100,1)), 44100)


%% ------------------------------- interpolate head-above-torso orientation

[l, r] = AKhrirInterpolation(0, 0, -50:1:50);

% quick plot (left ear only)
AKf
subplot(2,1,1)
    AKp(l, 't3d', 'y', -50:1:50, 'x', [0 3], 'dr', [-2 2])
    xlabel 'head-above-torso interpolation in degree'
subplot(2,1,2)
    AKp(l, 'm3d', 'y', -50:1:50, 'dr', [-20 20], 'N', 44100/100, 'x', [99 20000])
    xlabel 'head-above-torso interpolation in degree'
