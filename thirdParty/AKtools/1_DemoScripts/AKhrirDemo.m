% Two head-related transfer function (HRTF) databases are published by the
% Audio Communication Group:
% (1) The 'FABIAN HRTF database' contains HRTFs of the FABIAN head and
%     torso simulator for different head-above-torso orientations [1-3]
% (2) The 'HUTUBS HRTF database' contains HRTFs of 96 subjects [4-5]
%
% This script shows how to obtain head-related impulse responses (HRIRs)
% from the databases
%
%
% [1] Fabian Brinkmann, Alexander Lindau, Stefan Weinzierl, Steven van de
%     Par, Markus Müller-Trapet, Rob Opdam, and Michael Vorländer (2017) 
%     'A High Resolution and Full-Spherical Head-Related Transfer Function
%     Database for Different Head-Above-Torso Orientations', J. Audio Eng.
%     Soc, 65(10): 841--848. 
% [2] Fabian Brinkmann, Alexander Lindau, Stefan Weinzierl, Gunnar Geissler,
%     Steven van de Par, Markus Müller-Trapet, Rob Opdam, and Michael
%     Vorländer (2017) The FABIAN head-related transfer function data base. 
%     DOI: https://dx.doi.org/10.14279/depositonce-5718.2
% [3] F Brinkmann, R Roden,A Lindau, S Weinzierl "Audibility and inter-
%     polation of head-above-torso orientation in binaural technology."
%     IEEE J. Sel. Topics Signal Process., 9(5), 931-942, (2015).
% [4] tba
% [5] tba
%
% 2016/06 - fabian.brinkmann@tu-berlin.de
% 2018/19 - fabian.brinkmann@tu-berlin.de (added HUTUBS support)

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


%% ---------------------------------------- HUTUBS: check data availability
AKdependencies('HUTUBS')


%% ------------------------------------------------ HUTUBS: get single HRIR

% HRIR of participant 1 (FABIAN) for a source to the left
data = AKhrirDatabase('pp1_measured_sh', [90 0]);

% quick plot
AKf
AKp( squeeze(data), {'ir2d'; 'm2d'})
legend('left ear', 'right ear', 'location', 'SouthWest')


%% --------------------------------- HUTUBS: compare measured and simulated

% HRIR of participant 1 (FABIAN) for a source to the left
dataM = AKhrirDatabase('pp1_measured_sh',  [90 0]);
dataS = AKhrirDatabase('pp1_simulated_sh', [90 0]);

% quick plot
AKf
AKp( [dataM(:,1,1) dataS(:,1,1)], {'ir2d'; 'm2d'})
legend('left ear, measured', 'left ear, simulated', 'location', 'SouthWest')


%% ---------------------------------- HUTUBS: auralization of moving source

% diffuse field compensated auralization for participant 1 (FABIAN)
data = AKhrirDatabase('pp1_measured_auralization', [], true);

% uncomment the next line to listen via headphones
% soundsc(data, 44100);


%% ----------------------------------------- FABIAN: check data availablity

% The source position is interpolated using spherical harmonics
% coefficients. The HATO is interpolated using inverse distance weighting,
% as described in [3] using the magnitude and unwrapped phase spectra.

AKdependencies('FABIAN')


%% ---------------------------------------------------- FABIAN: single HRIR

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


%% -- FABIAN: single directional impuse response (diffuse field comensated)

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


%% ----------------------- FABIAN: interpolate head-above-torso orientation

[l, r] = AKhrirInterpolation(0, 0, -50:1:50);

% quick plot (left ear only)
AKf
subplot(2,1,1)
    AKp(l, 't3d', 'y', -50:1:50, 'x', [0 3], 'dr', [-2 2])
    xlabel 'head-above-torso interpolation in degree'
subplot(2,1,2)
    AKp(l, 'm3d', 'y', -50:1:50, 'dr', [-20 20], 'N', 44100/100, 'x', [99 20000])
    xlabel 'head-above-torso interpolation in degree'

    
%% --------------------------------------------- FABIAN vs. HUTUBS database

% The FABIAN head and torso simulator is contained in both databases
dataF = AKhrirInterpolation(90, 0, 0);
dataH = AKhrirDatabase('pp1_measured_ir', [90 0]);

% quick polot
AKf
AKp([dataF dataH(:,:,1)], {'ir2d'; 'm2d'})
legend('left ear, FABIAN database', 'left ear, HUTUBS database', 'location', 'SouthWest')