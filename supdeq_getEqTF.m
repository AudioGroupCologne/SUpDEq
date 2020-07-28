%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [eqTF_L, eqTF_R] = supdeq_getEqTF(eqDataset,samplingGrid,mode,channel,transformCore,phaseOnly)
%
% This function returns an equalization transfer function (separately for the left/right ear) 
% for an arbitrary direction, based on interpolation in the SH-domain and inverse 
% spherical Fourier transform.
%
% Output:
% eqTF_L / R    - Result of ISFT in Fourier domain (single-sided complex
%                 spectrum). Equalization transfer function (e.g., 
%                 rigid sphere transfer function) for L/R channel at each 
%                 spatial sampling point
%
% Input:
% eqDataset     - Struct with equalization dataset as SH-coefficients. 
%                 Can be the output of supdeq_getEqDataset.
% samplingGrid  - Spatial sampling grid (Q x 2 matrix), where the first 
%                 column holds the azimuth and the second the
%                 elevation (both in degree).
%                 Azimuth in degree (0=front, 90=left, 180=back, 270=right)
%                 (0 points to positive x-axis, 90 to positive y-axis)
%                 Elevations in degree (0=North Pole, 90=front, 180=South Pole)
%                 (0 points to positive z-axis, 180 to negative z-axis)
% mode          - 'DEG' or 'RAD' (radiant or degree).
%                 Default: 'DEG'
% channel       - Integer with 0, 1, or 2
%                 0 - Only left channel
%                 1 - Only right channel
%                 Default: 2 - Both channels (stereo)
% transformCore - String to define method to be used for the inverse 
%                 spherical Fourier transform. 
%                 'sofia - sofia_itc from SOFiA toolbox
%                 'ak'   - AKisht from AKtools 
%                 The results are exactly the same, but AKisht is faster 
%                 with big sampling grids
%                 Default: 'sofia'
% phaseOnly     - Set to 1 if eqTF_L/R should only contain the phase response 
%                 of eqDataset (allpass-fitlers) and not the magnitude response too
%                 Default: 0 (eqTF_L/R with magnitude and phase)
%
% Dependencies: SOFiA toolbox, AKtools
%
% Reference: -
%
% (C) 2020 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing 

function [eqTF_L, eqTF_R] = supdeq_getEqTF(eqDataset,samplingGrid,mode,channel,transformCore,phaseOnly)

if nargin < 3 || isempty(mode)
    mode = 'DEG';
end

if nargin < 4 || isempty(channel)
    channel = 2;
end

if nargin < 5 || isempty(transformCore)
    transformCore = 'sofia';
end

if nargin < 6 || isempty(phaseOnly)
    phaseOnly = 0;
end

%Get angles
az = samplingGrid(:,1);
el = samplingGrid(:,2);

%SOFiA needs RAD
if strcmp(transformCore,'sofia')
    if strcmp(mode,'DEG')
        az = az*pi/180;
        el = el*pi/180;
    end
end

%AK needs DEG
if strcmp(transformCore,'ak')
    if strcmp(mode,'RAD')
        az = az*180/pi;
        el = el*180/pi;
    end
end

%Check range
if sum(az < 0) > 0 || sum(el < 0) > 0
    error('Only positive azimuth/elevation values allowed!');
end

%% Perform transform

%Transform eqDataset to Fourier domain at sparse sampling grid points
%(inverse spherical Fourier transform)
[eqTF_L,eqTF_R] = supdeq_getArbHRTF(eqDataset,samplingGrid,mode,channel,transformCore);

if phaseOnly    
    %fprintf('Phase only eqTFs...\n')
    %Get only phase response of eqTF
    eqTF_L_phase = angle(eqTF_L);
    eqTF_R_phase = angle(eqTF_R);
    eqTF_L_mag = ones(size(eqTF_L,1),size(eqTF_L,2));
    eqTF_R_mag = ones(size(eqTF_R,1),size(eqTF_R,2));
    eqTF_L = eqTF_L_mag.*exp(1i*eqTF_L_phase);
    eqTF_R = eqTF_R_mag.*exp(1i*eqTF_R_phase);
end

end




