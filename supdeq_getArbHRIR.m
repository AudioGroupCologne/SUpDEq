%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [HRIR_L, HRIR_R] = supdeq_getArbHRIR(HRIRs_sfd,samplingGrid,mode,channel,transformCore)
%
% This function returns a HRIR / HRTF for an arbitrary direction, based on
% interpolation in the SH-domain and inverse spherical Fourier
% transform.
%
% Output:
% HRIR_L/R      - Result of ISFT and IFFT in time domain. HRIR for L/R
%                 channel at each spatial sampling point
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
%                  Default: 'DEG'
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
%
% Dependencies: SOFiA toolbox, AKtools
%
% Reference: 
% Bernsch�tz, B. (2016). Microphone Arrays and Sound Field Decomposition 
%                        for Dynamic Binaural Recording. TU Berlin. 
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH K�ln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing 

function [HRIR_L, HRIR_R] = supdeq_getArbHRIR(HRIRs_sfd,samplingGrid,mode,channel,transformCore)

if isfield(HRIRs_sfd,'FFToversize')
    FFToversize = HRIRs_sfd.FFToversize;
else
    FFToversize = 1;
end

if nargin < 3 || isempty(mode)
    mode = 'DEG';
end

if nargin < 4 || isempty(channel)
    channel = 2;
end

if nargin < 5 || isempty(transformCore)
    transformCore = 'sofia';
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

%% Perform transform with sofia_itc

if strcmp(transformCore,'sofia')
    
    if channel == 0 %Only left channel
        
        %Perform inverse spherical Fourier transform 
        Hl = sofia_itc(HRIRs_sfd.Hl_nm, [az el]);

        %Get mirror spectrum
        HlTwoSided = conj([Hl(:,:), conj(fliplr(Hl(:,2:end-1)))])';

        %Transform back to time domain
        HRIR_L = real(ifft(HlTwoSided));

        %Cut zeros depending on FFToversize
        finalLength = size(HRIR_L,1)/FFToversize;
        HRIR_L = HRIR_L(1:finalLength,:);
        HRIR_R = nan;
        
    elseif channel == 1 %Only right channel
        
        %Perform inverse spherical Fourier transform 
        Hr = sofia_itc(HRIRs_sfd.Hr_nm, [az el]);

        %Get mirror spectrum
        HrTwoSided = conj([Hr(:,:), conj(fliplr(Hr(:,2:end-1)))])';

        %Transform back to time domain
        HRIR_R = real(ifft(HrTwoSided));

        %Cut zeros depending on FFToversize
        finalLength = size(HRIR_R,1)/FFToversize;
        HRIR_L = nan;
        HRIR_R = HRIR_R(1:finalLength,:);
        
    elseif channel == 2 %Default - Stereo
        
        %Perform inverse spherical Fourier transform 
        Hl = sofia_itc(HRIRs_sfd.Hl_nm, [az el]);
        Hr = sofia_itc(HRIRs_sfd.Hr_nm, [az el]);

        %Get mirror spectrum
        HlTwoSided = conj([Hl(:,:), conj(fliplr(Hl(:,2:end-1)))])';
        HrTwoSided = conj([Hr(:,:), conj(fliplr(Hr(:,2:end-1)))])';

        %Transform back to time domain
        HRIR_L = real(ifft(HlTwoSided));
        HRIR_R = real(ifft(HrTwoSided));

        %Cut zeros depending on FFToversize
        finalLength = size(HRIR_L,1)/FFToversize;
        HRIR_L = HRIR_L(1:finalLength,:);
        HRIR_R = HRIR_R(1:finalLength,:);
        
    end
end

%% Perform transform with AKisht

if strcmp(transformCore,'ak')
    
    if channel == 0 %Only left channel
        
        %Perform inverse spherical Fourier transform 
        HRIR_L = AKisht(HRIRs_sfd.Hl_nm,true,[az el],'complex');

        %Cut zeros depending on FFToversize
        finalLength = size(HRIR_L,1)/FFToversize;
        HRIR_L = HRIR_L(1:finalLength,:);
        HRIR_R = nan;
        
    elseif channel == 1 %Only right channel
        
        %Perform inverse spherical Fourier transform 
        HRIR_R = AKisht(HRIRs_sfd.Hr_nm,true,[az el],'complex');

        %Cut zeros depending on FFToversize
        finalLength = size(HRIR_R,1)/FFToversize;
        HRIR_L = nan;
        HRIR_R = HRIR_R(1:finalLength,:);
        
    elseif channel == 2 %Default - Stereo

        %Perform inverse spherical Fourier transform 
        HRIR_L = AKisht(HRIRs_sfd.Hl_nm,true,[az el],'complex');
        HRIR_R = AKisht(HRIRs_sfd.Hr_nm,true,[az el],'complex');

        %Cut zeros depending on FFToversize
        finalLength = size(HRIR_L,1)/FFToversize;
        HRIR_L = HRIR_L(1:finalLength,:);
        HRIR_R = HRIR_R(1:finalLength,:);
        
    end
end

end




