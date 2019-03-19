%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [denseHRTFdataset, denseHRIRdataset, denseHRTFdataset_sfd] = supdeq_deq(eqHRTFdataset, deqDataset, N, samplingGrid, windowLength)
%
% This function performs the de-equalization of a equalized sparse HRTF 
% dataset, here called eqHRTFdataset, (in SH-domain / stored as SH-coefficients) 
% with the pre-defined deqDataset (dataset for de-equalization in SH-domain / 
% stored as SH-coefficients). The result is a dense hrtf dataset (according
% to the (dense) spatial sampling grid), saved as in structs as HRTFs
% (denseHRTFdataset), HRIRs (denseHRIRdataset) or as SH-coefficients
% (denseHRTFdataset_sfd)
%
% Output:
% denseHRTFdataset      - Dense HRTF dataset as a struct 
% denseHRIRdataset      - Dense HRIR dataset as a struct
% denseHRTFdataset_sfd  - Dense HRTF dataset as SH-coefficients,
%                         transformed with given order N
%
% Input:        
% eqHRTFdataset         - Struct with (SH-Coefficients) of the equalized 
%                         HRTF dataset for the left (Hl_nm) and right (Hr_nm) 
%                         channel/ear, absolute frequency scale f, 
%                         transform order N, and FFToversize
% deqDataset            - Struct with de-equalization dataset in SH-domain / as 
%                         SH-coefficients. Can be the output of supdeq_getEqDataset.
% N                     - Transform order N
% samplingGrid          - Spatial sampling grid (Q x 2 matrix), where the first 
%                         column holds the azimuth and the second the
%                         elevation (both in degree).
%                         Azimuth in degree
%                         (0=front, 90=left, 180=back, 270=right)
%                         (0 points to positive x-axis, 90 to positive y-axis)
%                         Elevations in degree
%                         (0=North Pole, 90=front, 180=South Pole)
%                         (0 points to positive z-axis, 180 to negative z-axis)
% windowLength          - Values of head and tail window in samples [a,b].
%                         Head and tail window will be a half sided hann
%                         window with the defined length. Windowing will be
%                         applied to the HRIRs. The HRTF dataset in frequency- 
%                         and SH-domain will be calculated based on these 
%                         windowed HRIRs. If the matrix is empty ([]), no
%                         windowing will be applied!
%                         Default: [], no windowing
%
% Dependencies: SOFiA toolbox
%
% (C) 2018 by CP,  Christoph Pörschmann
%             JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [denseHRTFdataset, denseHRIRdataset, denseHRTFdataset_sfd] = supdeq_deq(eqHRTFdataset, deqDataset, N, samplingGrid, windowLength)

if nargin < 5
    windowLength = [];
end

%Transform eqDataset to Fourier domain at dense sampling grid points
%(inverse spherical Fourier transform)
fprintf('Extracting %d deq transfer functions. This may take some time...\n',size(samplingGrid,1));
[deqTF_L,deqTF_R] = supdeq_getArbHRTF(deqDataset,samplingGrid);
fprintf('%d deq transfer functions extracted...\n',size(samplingGrid,1))

%Get HRTF from equalized HRTF dataset by inverse spherical Fourier transform at
%sampling points of dense sampling grid
fprintf('Extracting %d HRTFs. This may take some time...\n',size(samplingGrid,1));
[eqHRTF_L,eqHRTF_R] = supdeq_getArbHRTF(eqHRTFdataset,samplingGrid);
fprintf('%d HRTFs extracted...\n',size(samplingGrid,1))

%Perform de-equalization (multiplication in frequency domain)
HRTF_deequalized_L = zeros(size(samplingGrid,1),length(deqDataset.f));
HRTF_deequalized_R = zeros(size(samplingGrid,1),length(deqDataset.f));
for kk = 1:size(samplingGrid,1)
    HRTF_deequalized_L(kk,:)=eqHRTF_L(kk,:).*deqTF_L(kk,:);
    HRTF_deequalized_R(kk,:)=eqHRTF_R(kk,:).*deqTF_R(kk,:); 
    
    if ~mod(kk,100)
        fprintf('%d of %d HRTFs de-equalized...\n',kk,size(samplingGrid,1))
    end
    
end

%% Store in structs

% Transform HRTFs to time domain --> get HRIRs
%Get both sided spectrum
HRTF_deequalized_L_double = conj([HRTF_deequalized_L(:,:), conj(fliplr(HRTF_deequalized_L(:,2:end-1)))])';
HRTF_deequalized_R_double = conj([HRTF_deequalized_R(:,:), conj(fliplr(HRTF_deequalized_R(:,2:end-1)))])';
%Transform back to time domain
HRIR_deequalized_L = real(ifft(HRTF_deequalized_L_double));
HRIR_deequalized_R = real(ifft(HRTF_deequalized_R_double));
%Cut if FFToversize of eqHRTFdataset > 1
if eqHRTFdataset.FFToversize > 1
    newLength = size(HRIR_deequalized_L,1)/eqHRTFdataset.FFToversize;
    HRIR_deequalized_L = HRIR_deequalized_L(1:newLength,:);
    HRIR_deequalized_R = HRIR_deequalized_R(1:newLength,:);
end

%Window if windowLength is not empty
if ~isempty(windowLength) %windowing
    
    %Get window function in time domain
    win = supdeq_win(size(HRIR_deequalized_L,1),windowLength);
    %Apply to HRIRs
    for kk = 1:size(HRIR_deequalized_L,2)
        HRIR_deequalized_L(:,kk) = HRIR_deequalized_L(:,kk) .* win;
        HRIR_deequalized_R(:,kk) = HRIR_deequalized_R(:,kk) .* win;
    end
    
    %Store in struct denseHRIRdataset (window)
    denseHRIRdataset.HRIR_L = HRIR_deequalized_L'; %Flip again to have same array settings as in SOFiA
    denseHRIRdataset.HRIR_R = HRIR_deequalized_R';
    denseHRIRdataset.fs = eqHRTFdataset.f(end)*2;
    denseHRIRdataset.samplingGrid = samplingGrid;
    
    %Transform de-equalized windowed denseHRIRdataset to Fourier domain to get
    %denseHRTFdataset (with same FFToversize)
    HRTF_deequalized_win_L = fft(HRIR_deequalized_L',size(HRIR_deequalized_L,1)*eqHRTFdataset.FFToversize,2);
    HRTF_deequalized_win_L = HRTF_deequalized_win_L(:,1:end/2+1);
    HRTF_deequalized_win_R = fft(HRIR_deequalized_R',size(HRIR_deequalized_R,1)*eqHRTFdataset.FFToversize,2);
    HRTF_deequalized_win_R = HRTF_deequalized_win_R(:,1:end/2+1);
    
    %Store in struct denseHRTFdataset (window)
    denseHRTFdataset.HRTF_L = HRTF_deequalized_win_L;
    denseHRTFdataset.HRTF_R = HRTF_deequalized_win_R;
    denseHRTFdataset.f      = eqHRTFdataset.f;
    denseHRTFdataset.fs     = eqHRTFdataset.f(end)*2;
    denseHRTFdataset.FFToversize = eqHRTFdataset.FFToversize;
    denseHRTFdataset.samplingGrid = samplingGrid;
    
    % Transform de-equalized HRTFs to SH-domain --> get
    % denseHRTFdataset_sfd (window)
    fprintf('Performing spherical Fourier transform with N = %d. This may take some time...\n',N);
    denseHRTFdataset_sfd = supdeq_hrtf2sfd(HRTF_deequalized_win_L,HRTF_deequalized_win_R,N,samplingGrid);
    denseHRTFdataset_sfd.FFToversize = eqHRTFdataset.FFToversize;
    fprintf('Done with de-equalization...\n');
      
else %no windowing
    
    %Store in struct denseHRIRdataset (no window)
    denseHRIRdataset.HRIR_L = HRIR_deequalized_L'; %Flip again to have same array settings as in SOFiA
    denseHRIRdataset.HRIR_R = HRIR_deequalized_R';
    denseHRIRdataset.fs = eqHRTFdataset.f(end)*2;
    denseHRIRdataset.samplingGrid = samplingGrid;
    
    %Store in struct denseHRTFdataset (no window)
    denseHRTFdataset.HRTF_L = HRTF_deequalized_L;
    denseHRTFdataset.HRTF_R = HRTF_deequalized_R;
    denseHRTFdataset.f      = eqHRTFdataset.f;
    denseHRTFdataset.fs     = eqHRTFdataset.f(end)*2;
    denseHRTFdataset.FFToversize = eqHRTFdataset.FFToversize;
    denseHRTFdataset.samplingGrid = samplingGrid;

    % Transform de-equalized HRTFs to SH-domain --> get
    % denseHRTFdataset_sfd (no window)
    fprintf('Performing spherical Fourier transform with N = %d. This may take some time...\n',N);
    denseHRTFdataset_sfd = supdeq_hrtf2sfd(HRTF_deequalized_L,HRTF_deequalized_R,N,samplingGrid);
    denseHRTFdataset_sfd.FFToversize = eqHRTFdataset.FFToversize;
    fprintf('Done with de-equalization...\n');
    
end

end

