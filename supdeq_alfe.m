%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [eqHRTFDataset_alfe] = supdeq_alfe(eqHRTFDataset, samplingGrid, xover, fs)
%
% This function applies the ALFE (Adaptive Low Frequency Extention)
% algorithm in frequency domain (ALFE_fd) to HRTFs. Typically, it is 
% applied on the equalized HRTF dataset as a processing step of SUpDEq.
%
% Input:
% eqHRTFDataset         - Equalized HRTF dataset in the SH domain 
% xover                 - Crossover frequency in Hz
% fs                    - Sampling rate
%                         Default: 48000
% Output:       
% eqHRTFDataset_ALFE    - Equalized and ALFE-processeed HRTF dataset in the SH domain 
% 
% Dependencies: ALFE_fd
%
% References:
% Xie, Bosun: 'On the low frequency characteristics of head-related
% transfer function'. In: Chinese J. Acoust. 28(2):1-13 (2009).
%
% (C) 2019 by CP,  Christoph Pörschmann
%             JMA, Johannes M. Arend
%             Technische Hochschule Köln
%             University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [eqHRTFDataset_alfe] = supdeq_alfe(eqHRTFDataset, samplingGrid, xover, fs)

if nargin < 4 || isempty(fs)
    fs = 48000;
end

if length(xover)~=2
    warning('ALFE Crossover frequency wrong defined');
end

if or(min(xover)<300,max(xover)>1300) 
    warning('ALFE Crossover frequencies in not in typical range between 300 and 1300 Hz');
end

%Get equalized HRTFs from SH coefficients of equalized HRTF dataset
%[HRTF_equalized_L, HRTF_equalized_R] = supdeq_getArbHRTF(eqHRTFDataset,samplingGrid,'DEG',2,'ak');
%Get equalized HRTFs 
HRTF_equalized_L = eqHRTFDataset.HRTF_L;
HRTF_equalized_R = eqHRTFDataset.HRTF_R;
%Get mirror spectrum
HlTwoSided = AKsingle2bothSidedSpectrum(HRTF_equalized_L.');
HrTwoSided = AKsingle2bothSidedSpectrum(HRTF_equalized_R.');
%Transform back to time domain
hl = real(ifft(HlTwoSided));
hr = real(ifft(HrTwoSided));

%% Determine average value in crossover range

NFFT = size(hl,1);
n_link = round(xover / (fs/NFFT)) + 1;
ALFE_abs = mean(mean(abs(HlTwoSided(n_link(1):n_link(2),:)))+mean(abs(HrTwoSided(n_link(1):n_link(2),:))))/2;

%% Apply ALFE to all HRIRs

fprintf('Applying ALFE\n');

for kk = 1:length(hl(1,:))
    
    %Get ALFE processed HRIR
    hl_ALFE  = ALFE_fd(hl(:,kk),'f_link', xover,'fs',fs,'L_target',20*log10(ALFE_abs));
    hr_ALFE  = ALFE_fd(hr(:,kk),'f_link', xover,'fs',fs,'L_target',20*log10(ALFE_abs));

    %Interchange original HRIRs with ALFE processed HRIRs
    hl(:,kk)=hl_ALFE;
    hr(:,kk)=hr_ALFE;
    
    %Transform back to Fourier Domain
    Hl=fft(hl(:,kk));
    Hr=fft(hr(:,kk));
    
    %Get rid of mirror spectrum
    HRTF_equalized_ALFE_L(kk,:)=Hl(1:length(HRTF_equalized_L(1,:)));
    HRTF_equalized_ALFE_R(kk,:)=Hr(1:length(HRTF_equalized_R(1,:)));
    
    %Status prints
    if ~mod(kk,10)
        fprintf('|');
    end
    if ~mod(kk,500)
        fprintf(' %d\n',kk); 
    end
end

%Transform ALFE processed HRTFs to SH domain
fprintf('\nPerforming spherical Fourier transform with N = %d. This may take some time...\n',eqHRTFDataset.N);
eqHRTFDataset_temp = supdeq_hrtf2sfd(HRTF_equalized_ALFE_L,HRTF_equalized_ALFE_R,eqHRTFDataset.N,samplingGrid,fs,'ak');

%Write new output struct
eqHRTFDataset_alfe = eqHRTFDataset;
eqHRTFDataset_alfe.Hl_nm = eqHRTFDataset_temp.Hl_nm;
eqHRTFDataset_alfe.Hr_nm = eqHRTFDataset_temp.Hr_nm;
eqHRTFDataset_alfe.HRTF_L = HRTF_equalized_ALFE_L;
eqHRTFDataset_alfe.HRTF_R = HRTF_equalized_ALFE_R;
eqHRTFDataset_alfe.alfe = 1;
eqHRTFDataset_alfe.alfeXover = xover;
clear eqHRTFDataset_temp

fprintf('\nDone with ALFE...\n');
