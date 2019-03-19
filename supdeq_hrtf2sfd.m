%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function HRIRs_sfd = supdeq_hrtf2sfd(HRTF_L, HRTF_R, N, samplingGrid)
%
% This function transformes HRTFs to the SH-domain by spherical Fourier transform.
%
% Output:
% HRIRs_sfd     - Struct with Spherical Harmonic coefficients 
%                 (SH-Coefficients) for the left (Hl_nm) and right (Hr_nm) 
%                 channel/ear, absolute frequency scale f, 
%                 transform order N, and FFToversize
%
% Input:
% HRTF_L/R      - HRTF for left and right ear as [n x m] samples, with n =
%                 number of channels / sampling points and m = FFT bins
% N             - Transform order N
% samplingGrid  - Optional input [default = nan]
%                 If a spatial sampling grid is passed, it must be a Qx3 
%                 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
%                 In this case, the SOFiA functions "sofia_fdt" and 
%                 "sofia_stc" are used.
%                 If no samplingGrid (or []) is passed, the sampling grid
%                 according to the SOFA object is used for a spherical
%                 Fourier transform without sampling weights. In this case,
%                 the function AKsht is used.
%                 Azimuth positions in degree. Can be scalar or vector of size(el)
%                 (0=front, 90=left, 180=back, 270=right)
%                 (0 points to positive x-axis, 90 to positive y-axis)
%                 Elevations in degree. Can be scalar or vector of size(az)
%                 (0=North Pole, 90=front, 180=South Pole)
%                 (0 points to positive z-axis, 180 to negative z-axis)
% fs            - Sampling rate. Only needed to write frequency vector or
%                 if AKsht is used. 
%                 Default: 48000
%
% Dependencies: SOFiA toolbox, AKtools
%
% References:
% Benjamin Bernschütz: Microphone Arrays and Sound Field Decomposition 
% for Dynamic Binaural Recording. Ph.D. dissertation, Technical University
% Berlin (2016).
%
% Boaz Rafaely: Fundamentals of spherical array processing. In. Springer 
% topics in signal processing. Benesty, J.; Kellermann, W. (Eds.), 
% Springer, Heidelberg et al. (2015).
% 
% Franz Zotter: Analysis and synthesis of sound-radiation with spherical 
% arrays. Ph.D. dissertation, University of Music and Performing arts (2009).
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function HRIRs_sfd = supdeq_hrtf2sfd(HRTF_L, HRTF_R, N, samplingGrid, fs)

if nargin < 5
    fs = 48000;
end

%Check passed grid for weights  
weightsPassed = true;
if size(samplingGrid,2) ~= 3
    %warning('No sampling weights passed. Using AKsht to perform SHT with pseudo-inverse SH-matrix')
    weightsPassed = false;
end

%% If sampling grid with weights passed - use sofia_stc

if weightsPassed
        
    %Convert az and el to rad again
    samplingGrid(:,1:2) = samplingGrid(:,1:2) * pi / 180;
     
    %Transform to SH-domain with spherical Fourier transform (sofia_stc)
    %Get SH-coefficients for left and right channel
    Hl_nm = sofia_stc(N,HRTF_L,samplingGrid);
    Hr_nm = sofia_stc(N,HRTF_R,samplingGrid);

    %Write output struct
    HRIRs_sfd.Hl_nm         = Hl_nm;
    HRIRs_sfd.Hr_nm         = Hr_nm;
    HRIRs_sfd.f             = linspace(0,fs/2,size(HRTF_L,2));
    HRIRs_sfd.N             = N;
    HRIRs_sfd.FFToversize   = 1;
    
    disp('Transformation done with SOFiA toolbox');
end

%% If no weights passed - Use AKsht
if ~weightsPassed
    
    %Transform to SH-domain with spherical Fourier transform (AKsht)
    %Get SH-coefficients for left and right channel
    Hl_nm = AKsht(HRTF_L.',false,samplingGrid,N,'complex',fs);
    Hr_nm = AKsht(HRTF_R.',false,samplingGrid,N,'complex',fs);
    
    %Write output struct
    HRIRs_sfd.Hl_nm         = Hl_nm;
    HRIRs_sfd.Hr_nm         = Hr_nm;
    HRIRs_sfd.f             = linspace(0,fs/2,size(HRTF_L,2));
    HRIRs_sfd.N             = N;
    HRIRs_sfd.FFToversize   = 1;
    
    disp('Transformation done with AKTools');
end

end

