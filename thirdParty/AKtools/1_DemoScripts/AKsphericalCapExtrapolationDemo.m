% Demonstration of extrapolating missing data on a spherical sampling grid.
% In this case we assume that measured data is available above a given
% elevation, and all data below that elevation is extrapolated.
% The extrapolation is done using spherical harmonics as suggested by
% Ahrens et. al [1]:
%
% 1. A low order spherical harmonics transform is done using the available
%    data. The spherical harmonics representation of this data is then used
%    to calculate/extrapolate the missing data.
% 2. The available and missing data is joined to obtain full spherical data
% 3. A high order spherical harmonics transfrom is done using the joint
%    data to smooth the transition between the avaiable and extrapolated
%    data.
%
% [1] J. Ahrens, M.R.P. Thomas, I.J. Tashev (2012): "HRTF magnitude
%     modeling using a non-regularized least-squares fit of spherical
%     harmonics coefficients on incomplete data." APSIPA Annual Summit and
%     Conference, Hollywood, USA.
%
% 1/2015 - fabian.brinkmann@tu-berlin.de

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

% load HRIRs and spatial sampling grid
AKdependencies('FABIAN')

H            = SOFAload(fullfile('FABIAN_HRIR_modeled_HATO_0.sofa'));
hrir.raw     = shiftdim(H.Data.IR(:,1,:), 2);
grids.azFull = H.SourcePosition(:,1);
grids.elFull = H.SourcePosition(:,2);

clear H

%% ---------------------------------------- 1. set extrapolation parameters

N.high  = 35;          % High order
N.low   = 4;           % Low order
SHTmode = 'db_unwrap'; % SHT mode
elSplit = -64;         % Elevation below which the data will be interpolated


%% ------------------------------------- 2. split the spatial sampling grid
%  we assume that the data above elSplit is available and the data below is
%  missing and should be extrapolated
idA  = grids.elFull >= elSplit;
idM  = grids.elFull <  elSplit;

grids.az  = grids.azFull(idA);
grids.el  = grids.elFull(idA);
grids.azM = grids.azFull(idM);
grids.elM = grids.elFull(idM);

%% ----------------------------------- 3. get order limited reference HRIRs
%  we will use oder limited HRIRs - i.e. after spherical harmonics
%  transform, and inverse spherical harmonics transform - as a reference
%  for benchmarking the extrapolation. We use the full spharical data for
%  this purpose.

% precalculate SH matrixes
Ynm.high = AKsh(N.high, [], grids.azFull, 90-grids.elFull);
Ynm.low  = AKsh(N.low,  [], grids.azFull, 90-grids.elFull);

% get HRIRs
[fnm.ref, ~, isEven] = AKsht(hrir.raw, true, pinv(Ynm.high), N.high, SHTmode);
hrir.ref = AKisht(fnm.ref, true, Ynm.high, SHTmode, isEven);


%% ------------------------------------------------------- 4. extrapolation

% Apply a low order spherical harmonics transfrom to the available data.
fnm.low     = AKsht(hrir.ref(:,idA), true, pinv(Ynm.low(idA,:)), N.low, SHTmode);
hrir.low    = AKisht(fnm.low, true, Ynm.low, SHTmode, isEven);
condYnm.low = cond(Ynm.low(idA,:));

% Use the low order fit to complete the missing data
hrir.cat = [hrir.raw(:,idA) hrir.low(:,idM)];

% get HRIRs including spherical cap extrapolation from high order spherical
% harmonics transfrom (This smoothes the transition between available and
% extrapolated data)
fnm.high     = AKsht(hrir.cat, true, pinv(Ynm.high), N.high, SHTmode);
hrir.high    = AKisht(fnm.high, true, Ynm.high, SHTmode, isEven);
condYnm.high = cond(Ynm.high);

% get difference to original
hrir.diff = ifft(abs(fft(hrir.ref)) ./ abs(fft(hrir.high)), 'symmetric');


%% ---------------------------------------------------------------- 4. plot
%  we look at the magnitude spectra in the sagittal median plane

% get median plane
id = AKsubGrid([[grids.az; grids.azM] [grids.el; grids.elM]], 'sagittal', 0);
id = circshift(id, numel(id)/2);

AKf(40, 10)
subplot(1,5,1)
    AKp(hrir.ref(:,id), 'm3d', 'dr', [-20 20], 'cb', 'south', 'cm', 'AKbwr')
    title('a: Reference data')
subplot(1,5,2)
     AKp(hrir.low(:,id), 'm3d', 'dr', [-20 20], 'cb', 'south', 'cm', 'AKbwr')
     title('b: Low order fit')
subplot(1,5,3)
     AKp(hrir.cat(:,id), 'm3d', 'dr', [-20 20], 'cb', 'south', 'cm', 'AKbwr')
     title('c: ''a'' and ''b'' joined')
subplot(1,5,4)
     AKp(hrir.high(:,id), 'm3d', 'dr', [-20 20], 'cb', 'south', 'cm', 'AKbwr')
     title('d: ''c'' after SHT, and ISHT')
subplot(1,5,5)
     AKp(hrir.diff(:,id), 'm3d', 'dr', [-10 10], 'cb', 'south', 'cm', 'AKbwr')
     title('e: difference between ''a'', and ''d''')

     
%% -------------------------------------------------------- 5. auralization
%  auralization of the median plane HRIRs. We will assume that the left and
%  right ear HRIRs are identical (we only extrapolated the data for the
%  right ear)

% render the auralization using pink noise
y.ref  = AKshAura(AKnoise(5*44100,1), [0 0 NaN 180 180], [-90 90 NaN 90 -90], fnm.ref,  fnm.ref);
y.low  = AKshAura(AKnoise(5*44100,1), [0 0 NaN 180 180], [-90 90 NaN 90 -90], fnm.low,  fnm.low);
y.high = AKshAura(AKnoise(5*44100,1), [0 0 NaN 180 180], [-90 90 NaN 90 -90], fnm.high, fnm.high);

% listen

% soundsc(y.ref,  44100)
% soundsc(y.low,  44100)
% soundsc(y.high, 44100)
