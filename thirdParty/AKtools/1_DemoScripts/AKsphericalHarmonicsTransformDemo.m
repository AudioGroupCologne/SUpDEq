% Demonstrates basic usage of the spherical harmonics transform. Function
% definition follows [1]. The transformation is done using the pseudo
% inverse of the matrix containing the SH basis functions at the positions
% of the sampling grid.
%
% [1] Boaz Rafaely: Fundamentals of spherical array processing. In.
% Springer topics in signal processing. Benesty, J.; Kellermann, W. (Eds.),
% Springer, Heidelberg et al., first edition (2015).
%
% 09/2015    - fabian.brinkmann@tu-berlin.de,

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

%#ok<*NASGU>
close all; clear; clc

% suppress unused variables warning for demo script
%#ok<*ASGLU>

% get data for testing spherical harmonics transform
% (we use HRIRs for this)
AKdependencies('FABIAN')

H  = SOFAload('FABIAN_HRIR_measured_HATO_0.sofa');
fs = H.Data.SamplingRate;
g  = H.SourcePosition(:, 1:2);
H  = shiftdim(H.Data.IR, 2);
l  = H(:,:,1);
r  = H(:,:,2);

clear H


% definition of elevation is different for spherical harmonics:
% 0 degree: north pole, 90 degree: front, 180 degree south pole
g(:,2) = 90-g(:,2);

%% ------------------------------------------------- 1. transform parameter

% SH order
% This determines the spatial resolution of the SH transform (The higher
% the order, the higher the resolution). If the order is to high, you will
% obtain spatial aliasing. If it is to low you will have significant
% truncation error, i.e. loss of data due to the limited spatial
% resolution.
% With this data you shoud be able to go up to N = 35, we start with N = 10
% to save time...
N       = 10;

% Select the input data for the SH transform
% 'complex'    - transform on complex spectra
% 'abs_unwrap' - transform on magnitude, and unwrapped phase spectrum
% 'db_unwrap'  - transform on log. magnitude, and unwrapped phase spectrum
%
% Note: 'abs_unwrap', and 'db_unwrap' are prone to errors in phase
% unwrapping. This might occur if the phase at low frequencies is noisy.
% However, this yields better results for compared to 'complex', at least
% for low truncation orders and the contra lateral ear
SHTmode = 'complex';

% Handling of SHT mode
% compact = false - separate transfrom on magnitude and phase spectra. This
%                   means that you will obtain two sets of SH coefficients
%                   (one for magnitude, and one for phase)
% compact = true  - transfrom on magnitude + 1j unwrapped phase. This means
%                   you obtain only one set of coefficients
%
% Note that this does not change the results after inverse transfrom, but
% in some cases you might want to manipulate magnitude and phase data
% seperately.
compact = true;

% get matrices with SH basis function
Ynm = AKsh(N, [], g(:,1), g(:,2));

% check the condition of the matrix - it should be close to 1
YnmCond = cond(Ynm);
disp(['SH matrix condition: ' num2str(round(YnmCond*100)/100)])

% get the inverse for the SH transform
YnmInv = pinv(Ynm);


%% ---------------------------------- 2: SH transform and inverse transfomr

% SH transform
[f_nm.l, f, is_even] = AKsht(l, true, YnmInv, N, SHTmode, fs, compact);
f_nm.r               = AKsht(r, true, YnmInv, N, SHTmode, fs, compact);
% Note: You can also run the transfrom without using the inverted Matrix. If
% you do this, it is calculated inside AKsht.m wich takes longer if you run
% multiple transforms on the identical grid:
% [f_nm.l, f, is_even] = AKsht(l, true, g, N, SHTmode, fs, compact);

% inverse SH transform
li = AKisht(f_nm.l, true, Ynm, SHTmode, is_even, compact);
ri = AKisht(f_nm.r, true, Ynm, SHTmode, is_even, compact);
% Note: You can again pass the grid instead of the Ynm matrix
% li = AKisht(f_nm.l, true, g, SHTmode, is_even, compact);


%% ---------------------------------- 3. plot before and after SH transform

% choose the source position that you want to see [azimuth elevation]
plotAngle = [90 0];  % <- this is a source to the left on the horizontal plane
id = AKsubGrid([g(:,1) 90-g(:,2)], 'any', plotAngle);

% plot
AKf
subplot(2,2,1)
    AKp(l(:,id), 't2d')
    AKp(li(:,id), 't2d', 'c', 'r')
    title 'Left ear HRIR'
    legend('original', 'sh transformed', 'location', 'NorthEast')
subplot(2,2,2)
    AKp(r(:,id), 't2d')
    AKp(ri(:,id), 't2d', 'c', 'r')
    title 'Right ear HRIR'
subplot(2,2,3)
    AKp(l(:,id), 'm2d')
    AKp(li(:,id), 'm2d', 'c', 'r')
    title 'Left ear HRTF'
subplot(2,2,4)
    AKp(r(:,id), 'm2d')
    AKp(ri(:,id), 'm2d', 'c', 'r')
    title 'Right ear HRTF'

%% ------------------------------ 4. auralization (watch your audio level!)

% moving source on the horizontal plane
y = AKshAura(AKnoise(5*fs), [0 360], [0 0], f_nm.l, f_nm.r, SHTmode, is_even);

% soundsc(y*.5, fs);

%% moving source on the median plane
y = AKshAura(AKnoise(5*fs), [0 0 NaN 180 180], [90 -90 NaN -90 90], f_nm.l, f_nm.r, SHTmode, is_even);

% soundsc(y*.5, fs);

%% moving source on the frontal plane
y = AKshAura(AKnoise(5*fs), [90 90 NaN 270 270], [90 -90 NaN -90 90], f_nm.l, f_nm.r, SHTmode, is_even);

% soundsc(y*.5, fs);
