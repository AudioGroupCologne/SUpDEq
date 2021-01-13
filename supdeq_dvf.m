%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [shiftedHRTF_L, shiftedHRTF_R, ffShiftedHRTFdataset] = supdeq_dvf(HRTFdataset, newDistance, samplingGrid ,shiftToFarField, adjustAmplitude, radius, earPosition, N_HRTF)
%
% This function applies a distance variation shift to HRTFs using distance 
% variation functions (DVF). As it works with HRTF in SH domain, acoustic 
% parallax effects (high frequency parallax effects induced by the pinna 
% which are not covered by the rigid sphere transfer functinos) can be
% embedded into the distance shifted HRTFs.
%
% Output:
% shiftedHRTF_L/R      - Shifted HRTF_L/R in frequency domain (single sided complex spectrum)
% ffShiftedHRTFdataset - Struct with far-field shifted HRTF dataset in SH
%                        domain according to input HRTFdataset. Only if
%                        'shiftToFarField' was applied
%
% Input:
% HRTFdataset           - Struct with the HRTF dataset as
%                         SH-coefficients for the left (Hl_nm) and 
%                         right (Hr_nm) channel/ear, absolute frequency scale f, 
%                         transform order N, FFToversize and sourceDistance
%                         If field 'sourceDistance' does not exist, the
%                         function assumes far-field measurements and option
%                         'shiftToFarField' will always be set to 'false'
% newDistance           - New distance of HRTF dataset in m 
% samplingGrid          - Spatial sampling grid (Q x 2 matrix) for the distance shifted HRTFs, 
%                         where the first column holds the azimuth and the second 
%                         the elevation (both in degree).
%                         Azimuth in degree (0=front, 90=left, 180=back, 270=right)
%                         (0 points to positive x-axis, 90 to positive y-axis)
%                         Elevations in degree (0=North Pole, 90=front, 180=South Pole)
%                         (0 points to positive z-axis, 180 to negative z-axis)
% shiftToFarField       - Boolean with true / false. Option to shift
%                         measurements, which strictly speaking were not measured in 
%                         the far-field (e.g., at d = 3.25 m as our Neumann KU100 set) 
%                         and thus still have acoustic parallax effects, to far field before 
%                         shifting to new distance. 
%                         Default: false
%                         Note: If set to false, only the inverse SH
%                         transform will be applied to get HRTFs for the
%                         new distance according to the sampling grid.
%                         Thus, processing will be much faster. If set to
%                         true, various forward and inverse SH transform
%                         are gonne be performed on a high order grid, so
%                         processing takes longer.
% adjustAmplitude       - Boolen with true / false. Option to adjust the HRTF level
%                         according to the 1/r law.
%                         Default: false --> HRTF have roughly the same
%                         amplitude after distance shift
% radius                - Define radius of rigid sphere / HRTFs in m
%                         Default: 0.0875
% earPosition           - 4 x 1 row vector describing the position of the ears in
%                         spherical coordinates in degree [azL, elL, azR, elR]
%                         Default: [90, 90, 270, 90] (left-right symmetrical)
% N_HRTF                - Spherical harmonics order for HRTF processing. If
%                         not specified, Nmax from the HRTFdataset struct will be applied (mostly
%                         the best choice)
%                         Default: Nmax from HRTFdataset struct
%
% Dependencies: -
%
% References:
% [1] A. Kan, C. Jin, and A. van Schaik, “A psychophysical evaluation of 
% near-field head-related transfer functions synthesized using a distance 
% variation function,” J. Acoust. Soc. Am., vol. 125, no. 4, pp. 2233–2242, 2009.
%
% [2] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and 
% perception of nearby virtual sound sources with a near-field filter model,” 
% Appl. Acoust., vol. 115, pp. 61–73, 2017.
%
% [3] D. Romblom and B. Cook, “Near-Field Compensation for HRTF Processing,” 
% in Proceedings of the 125th AES Convention, San Francisco, USA, 2008, pp. 1–6.
%
% (C) 2020 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [shiftedHRTF_L, shiftedHRTF_R, ffShiftedHRTFdataset] = supdeq_dvf(HRTFdataset, newDistance, samplingGrid, shiftToFarField, adjustAmplitude, radius, earPosition, N_HRTF)

if nargin < 2 || isempty(newDistance)
   error('Please specify new distance in m'); 
end

if nargin < 3 || isempty(samplingGrid)
   error('Please define sampling grid'); 
end

if nargin < 4 || isempty(shiftToFarField)
   shiftToFarField = false;
end

if nargin < 5 || isempty(adjustAmplitude)
   adjustAmplitude = false; 
end

if nargin < 6 || isempty(radius)
   radius = 0.0875; 
end

if nargin < 7 || isempty(earPosition)
    earPosition = [90, 90, 270, 90];
end

if nargin < 8 || isempty(N_HRTF)
    if isfield(HRTFdataset,'N')
        N_HRTF = HRTFdataset.N; 
    end
    if isfield(HRTFdataset,'Nmax')
        N_HRTF = HRTFdataset.Nmax; 
    end
end

if isfield(HRTFdataset,'sourceDistance')
    md = HRTFdataset.sourceDistance;
else
    md = 100; %Assume far field
    shiftToFarField = false;
end

NFFT = length(HRTFdataset.f)*2-2;
fs = HRTFdataset.f(end)*2;
ffShiftedHRTFdataset = nan;

%% Get various grids

%Omit weights of samplingGrid if defined
samplingGrid = samplingGrid(:,1:2);

%Get sg for STF
sgSTF = samplingGrid; sgSTF(:,2) = 90-sgSTF(:,2);
sgSTF_nd = sgSTF; sgSTF_nd(:,3) = newDistance;
sgSTF_md = sgSTF; sgSTF_md(:,3) = md;
ffDistance = 100;
sgSTF_ff = sgSTF; sgSTF_ff(:,3) = ffDistance;

earPositionSTF = earPosition; earPositionSTF(2) = 90-earPositionSTF(2); earPositionSTF(4) = 90-earPositionSTF(4);

%Get parallax grids of new distance
[sg_nd_para_l,sg_nd_para_r] = supdeq_parallax(samplingGrid,newDistance,radius,earPosition);

%% Get rigid sphere transfer functions for new distance

stf_nd      = AKsphericalHead(sgSTF_nd,earPositionSTF,false,radius,newDistance,100,NFFT,fs);
stf_nd_L    = AKboth2singleSidedSpectrum(fft(squeeze(stf_nd(:,:,1))));
stf_nd_R    = AKboth2singleSidedSpectrum(fft(squeeze(stf_nd(:,:,2))));


%% If dataset should be shifted to far field first before shifting to new distance
%This compensates parallax effects which are already in the measurement

if shiftToFarField
    
    disp('Shifting measured HRTFs first to far field');
    
    %Get Lebedev sampling grid according to N
    lebN = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65];
     
    %Get closest N
    [~,cId] = min(abs(N_HRTF-lebN));
    N_HRTF_leb = lebN(cId);
        
    if N_HRTF_leb < N_HRTF
        cId = cId + 1;
        N_HRTF = lebN(cId);
    else
        N_HRTF = N_HRTF_leb;
    end

    %Get sampling grid according to N_HRTF for internal sh processing
    sg_int = supdeq_lebedev([],N_HRTF); 
    %Omit weights
    sg_int = sg_int(:,1:2);
    
    %Get parallax grid for measurement distance
    [sg_md_para_l,sg_md_para_r] = supdeq_parallax(sg_int,md,radius,earPosition);
    
    %Get STFs for far field position on dense grid
    ffDistance = 100;
    sgSTF_int_ff = sg_int; sgSTF_int_ff(:,2) = 90-sgSTF_int_ff(:,2);
    sgSTF_int_ff(:,3) = ffDistance;
    stf_int_ff  = AKsphericalHead(sgSTF_int_ff,earPositionSTF,false,radius,ffDistance,100,NFFT,fs);
    stf_int_ff_L= AKboth2singleSidedSpectrum(fft(squeeze(stf_int_ff(:,:,1))));
    stf_int_ff_R= AKboth2singleSidedSpectrum(fft(squeeze(stf_int_ff(:,:,2))));
    
    %Get STFs for measurement position on dense grid
    sgSTF_int_md = sg_int; sgSTF_int_md(:,2) = 90-sgSTF_int_md(:,2);
    sgSTF_int_md(:,3) = md;
    stf_int_md  = AKsphericalHead(sgSTF_int_md,earPositionSTF,false,radius,md,100,NFFT,fs);
    stf_int_md_L= AKboth2singleSidedSpectrum(fft(squeeze(stf_int_md(:,:,1))));
    stf_int_md_R= AKboth2singleSidedSpectrum(fft(squeeze(stf_int_md(:,:,2))));
    
    %Calculate DVFs for left and right ear, Target = far
    DVF_int_L = stf_int_ff_L./stf_int_md_L; 
    DVF_int_R = stf_int_ff_R./stf_int_md_R; 
    
    %Get parallaxe grid HRTFs of measured dataset for measurement position
    [hrtf_L] = supdeq_getArbHRTF(HRTFdataset,sg_md_para_l,'DEG',0,'ak');
    [~,hrtf_R] = supdeq_getArbHRTF(HRTFdataset,sg_md_para_r,'DEG',1,'ak');
    
    %Apply DVF to shift to far-field
    ff_hrtf_L = (DVF_int_L.').*hrtf_L; 
    ff_hrtf_R = (DVF_int_R.').*hrtf_R; 
    
    %Transform far-field shifted HRTF set to SH domain again for further
    %processing
    ffShiftedHRTFdataset = supdeq_hrtf2sfd(ff_hrtf_L,ff_hrtf_R,N_HRTF,sg_int,fs,'ak');
    
    %Get some parameters
    ffShiftedHRTFdataset.FFToversize = HRTFdataset.FFToversize;
    ffShiftedHRTFdataset.sourceDistance = 100;
    ffShiftedHRTFdataset.ffShifted = true;
end

%% Apply DVFs

disp('Shifting HRTFs to new position');

if shiftToFarField
    %Calculate DVFs for left and right ear, Target = close
    %Work with far-field distance STFs
    stf_ff      = AKsphericalHead(sgSTF_ff,earPositionSTF,false,radius,ffDistance,100,NFFT,fs);
    stf_ff_L    = AKboth2singleSidedSpectrum(fft(squeeze(stf_ff(:,:,1))));
    stf_ff_R    = AKboth2singleSidedSpectrum(fft(squeeze(stf_ff(:,:,2))));
    
    DVF_L = stf_nd_L./stf_ff_L; 
    DVF_R = stf_nd_R./stf_ff_R; 
    
    %Get parallaxe grid HRTFs of far-field dataset for new position 
    [hrtf_L] = supdeq_getArbHRTF(ffShiftedHRTFdataset,sg_nd_para_l,'DEG',0,'ak');
    [~,hrtf_R] = supdeq_getArbHRTF(ffShiftedHRTFdataset,sg_nd_para_r,'DEG',1,'ak');
    
else
    %Calculate DVFs for left and right ear, Target = close
    %Work with measurement distance STFs
    stf_md      = AKsphericalHead(sgSTF_md,earPositionSTF,false,radius,md,100,NFFT,fs);
    stf_md_L    = AKboth2singleSidedSpectrum(fft(squeeze(stf_md(:,:,1))));
    stf_md_R    = AKboth2singleSidedSpectrum(fft(squeeze(stf_md(:,:,2))));
    
    DVF_L = stf_nd_L./stf_md_L; 
    DVF_R = stf_nd_R./stf_md_R; 
    
    %Get parallaxe grid HRTFs of far-field dataset for new position 
    [hrtf_L] = supdeq_getArbHRTF(HRTFdataset,sg_nd_para_l,'DEG',0,'ak');
    [~,hrtf_R] = supdeq_getArbHRTF(HRTFdataset,sg_nd_para_r,'DEG',1,'ak');
end

%Apply DVFs
shiftedHRTF_L = (DVF_L.').*hrtf_L; 
shiftedHRTF_R = (DVF_R.').*hrtf_R; 

if adjustAmplitude
    
    %Set referenceDistance to 1 if md was not defined
    if md == 100
        refDistance = 1;
    else
        refDistance = md;
    end
    
    shiftedHRTF_L = shiftedHRTF_L*(refDistance/newDistance);
    shiftedHRTF_R = shiftedHRTF_R*(refDistance/newDistance);
    
end

disp('Done...');

end

