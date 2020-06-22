% Demonstrates how to obtain head-related impulse responses (HRIR) of
% the spherical head model.
%
% Additional information on the model itself can be found in
% 2_Tools/SphericalHarmonics/AKsphericalHead.pdf
%
% 12/2017 - fabian.brinkmann@tu-berlin.de

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

%% -------------------------- 1. HRIR for a source position at the left ear
%  we only specify the source position, and use the default parameter for
%  everything else

h = AKsphericalHead([90 0]);

AKf(20,20)
AKp( squeeze( h ), {'t2d';'m2d'} )
legend('left ear', 'right ear', 'location', 'SouthWest')

%% ------------------------------------------ 2. HRIRs at varying distances

h = AKsphericalHead([90 0 2; 90 0 1; 90 0 .5]);

AKf(20,10)
subplot(1,2,1)
    AKp(h(:,:,1), 'm2d', 'dr', [-20 20])
    title 'left ear'
subplot(1,2,2)
    AKp(h(:,:,2), 'm2d', 'dr', [-20 20])
    title 'right ear'
    legend('2 m', '1 m', '0.5 m', 'location', 'SouthWest')

    
%% --------------------------------------- 3. HRIRs in the horizontal plane

[~, sg] = AKsubGrid(2, 'transverse', 0);

h = AKsphericalHead(sg);

AKf(30,20)
AKp(h(:,:,1), {'et3d', ''}, 'dr', [10 -40], 'y', sg(:,1), 'x', [0 5])
title 'left ear HRIRs (log)'; xlabel 'azimuth in deg.'
AKp(h(:,:,1), {'', 'm3d'}, 'dr', [10 -20], 'y', sg(:,1))
title 'left ear HRTF'; 
xlabel 'azimuth in deg.'
    

%% -------------- 4. HRIRs in the median plane with and without offset ears

[~, sg, sgPol] = AKsubGrid(2, 'sagittal', 0);

h = AKsphericalHead(sg);
h_symmetric = AKsphericalHead(sg, [90 0]);

AKf(30,20)
AKp(h(:,:,1), {'et3d' ''; '' ''}, 'y', sgPol(:,2), 'dr', [10 -40], 'x', [0 5])
title 'left ear HRIRs (offset ears)'; xlabel 'elevation in deg.'
AKp(h(:,:,1), {'' ''; 'm3d' ''}, 'y', sgPol(:,2), 'dr', [10 -10])
title 'left ear HRTFs (offset ears)'; xlabel 'elevation in deg.'

AKp(h_symmetric(:,:,1), {'' 'et3d'; '' ''}, 'y', sgPol(:,2), 'dr', [10 -40], 'x', [0 5])
title 'left ear HRIRs (symmetric ears)'; xlabel 'elevation in deg.'
AKp(h_symmetric(:,:,1), {'' ''; '' 'm3d'}, 'y', sgPol(:,2), 'dr', [10 -10])
title 'left ear HRTFs (symmetric ears)'; xlabel 'elevation in deg.'

%% ------------------------------------------ 5. Diffuse field compensation
%  in case the spherical head model is used for a headphone based
%  auralization, it might be a good idea to compensate for the diffuse
%  field transfer function of the spherical head

% spherical head directional IR
h = AKsphericalHead([90 0]);

% spherical head diffuse IR
[hDiffuse, hDiffuseInverted] = AKsphericalHeadDiffuse;

% apply the diffuse field compensation
hCompensated = fftfilt(hDiffuseInverted, h(:,:,1));

AKf(20,10)
AKp([h(:,:,1) hCompensated], 'm2d')
AKp(hDiffuseInverted, 'm2d', 'N', 1024)
title 'magnitude spectrum (left ear)'
legend('original', 'compensated', 'compensation filter', 'location', 'SouthWest')