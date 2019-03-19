%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function hrtf = supdeq_getArbHRTF(HRIRs_sfd,samplingGrid,mode)
%
% This function returns a HRTF (separately for the left/right ear) for an 
% arbitrary direction, based on interpolation in the SH-domain and inverse 
% spherical Fourier transform.
%
% Output:
% HRTF_L / R    - Result of ISFT in Fourier domain (single-sided complex
%                 spectrum). HRTF for L/R channel at each spatial sampling point
%
% Input:
% HRIRs_sfd     - Struct with HRIRs in spatial Fourier domain (sfd) / SH-domain
% samplingGrid  - Spatial sampling grid (Q x 2 matrix), where the first 
%                 column holds the azimuth and the second the
%                 elevation (both in degree).
%                 Azimuth in degree (0=front, 90=left, 180=back, 270=right)
%                 (0 points to positive x-axis, 90 to positive y-axis)
%                 Elevations in degree (0=North Pole, 90=front, 180=South Pole)
%                 (0 points to positive z-axis, 180 to negative z-axis)
% mode          - 'DEG' or 'RAD' (radiant or degree).
%                 Default: 'DEG'
%
% Dependencies: SOFiA toolbox
%
% Reference: 
% Bernschütz, B. (2016). Microphone Arrays and Sound Field Decomposition 
%                        for Dynamic Binaural Recording. TU Berlin.
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing 

function [HRTF_L, HRTF_R] = supdeq_getArbHRTF(HRIRs_sfd,samplingGrid,mode)

if nargin < 3
    mode = 'DEG';
end

%Get angles
az = samplingGrid(:,1);
el = samplingGrid(:,2);

%SOFiA needs RAD
if strcmp(mode,'DEG')
    az = az*pi/180;
    el = el*pi/180;
end

%Check range
if sum(az < 0) > 0 || sum(el < 0) > 0
    error('Only positive azimuth/elevation values allowed!');
end

%%
%Perform inverse spherical Fourier transform 
HRTF_L = sofia_itc(HRIRs_sfd.Hl_nm, [az el]);
HRTF_R = sofia_itc(HRIRs_sfd.Hr_nm, [az el]);

end




