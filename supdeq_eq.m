%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [eqHRTFdataset, HRTF_equalized_L, HRTF_equalized_R] = supdeq_eq(sparseHRTFdataset, eqDataset, N, samplingGrid, tikhEps)
%
% This function performs the equalization of a sparse HRTF dataset (in
% SH-domain / stored as SH-coefficients) with the pre-defined eqDataset
% (dataset for equalization in SH-domain / stored as SH-coefficients). The
% result, the eqHRTFdataset, is the equalized counterpart of the input HRTF
% dataset, stored as SH-coefficients. Transform order N of the
% eqHRTFdataset is set to the maximum order with respect to the sparse
% sampling grid. For example, a lebedev grid with 86 nodes has Nmax = 7.
%
% Output:
% eqHRTFdataset         - Struct with the equalized HRTF dataset as
%                         SH-coefficients for the left (Hl_nm) and 
%                         right (Hr_nm) channel/ear, equalized HRTFs
%                         for the left (HRTF_L) and right (HRTF_R) ear,
%                         absolute frequency scale f, transform order N, 
%                         and FFToversize
% HRTF_equalized_L/R    - Equalized complex HRTFs 
%
% Input:        
% sparseHRTFdataset     - Struct with sparse HRTF dataset in frequency
%                         domain ('normal' HRTFs...). Struct needs to
%                         provide HRTF_L/HRTF_R, the samplingGrid, and Nmax
% eqDataset             - Struct with equalization dataset as 
%                         SH-coefficients. Can be the output of 
%                         supdeq_getEqDataset.
% N                     - Transform order N
% samplingGrid          - Spatial sampling grid (Q x 2 or Q x 3 matrix), 
%                         where the first column holds the azimuth, 
%                         the second the elevation (both in degree), and
%                         optionally the third the sampling weights.
%                         Azimuth in degree
%                         (0=front, 90=left, 180=back, 270=right)
%                         (0 points to positive x-axis, 90 to positive y-axis)
%                         Elevations in degree
%                         (0=North Pole, 90=front, 180=South Pole)
%                         (0 points to positive z-axis, 180 to negative z-axis)
% tikhEps               - Define epsilon of Tikhonov regularization if
%                         regularization should be applied
%                         Applying the Tikhonov regularization will always result 
%                         in a least-square fit solution for the SH transform. 
%                         Variable 'transformCore' is neglected when 'tikhEps' is 
%                         defined as the regularized least-square spherical Fourier 
%                         transform is applied directly without any third party 
%                         toolbox. Depending on the sampling grids, weights are
%                         applied or not.
%                         Default: 0 (no Tikhonov regularization)
%
% Dependencies: SOFiA toolbox, AKTools
%   
% (C) 2018/2019 by  JMA, Johannes M. Arend
%                   CP,  Christoph Pörschmann
%                   TH Köln - University of Applied Sciences
%                   Institute of Communications Engineering
%                   Department of Acoustics and Audio Signal Processing

function [eqHRTFdataset, HRTF_equalized_L, HRTF_equalized_R] = supdeq_eq(sparseHRTFdataset, eqDataset, N, samplingGrid, tikhEps)

if nargin < 5 || isempty(tikhEps)
    tikhEps = 0;
end

%Get fs
fs = sparseHRTFdataset.f(end)*2;

%Transform eqDataset to Fourier domain at sparse sampling grid points
%(inverse spherical Fourier transform)
fprintf('Extracting %d eq transfer functions. This may take some time...\n',size(samplingGrid,1));
[eqTF_L,eqTF_R] = supdeq_getArbHRTF(eqDataset,samplingGrid,[],[],'ak'); %Use AKtools...
fprintf('%d eq transfer functions extracted...\n',size(samplingGrid,1))

%Perform equalization (division in frequency domain)
HRTF_equalized_L = zeros(size(samplingGrid,1),length(eqDataset.f));
HRTF_equalized_R = zeros(size(samplingGrid,1),length(eqDataset.f));
for kk = 1:size(samplingGrid,1)
    HRTF_equalized_L(kk,:)=sparseHRTFdataset.HRTF_L(kk,:)./eqTF_L(kk,:);
    HRTF_equalized_R(kk,:)=sparseHRTFdataset.HRTF_R(kk,:)./eqTF_R(kk,:); 
    
    if ~mod(kk,10)
        fprintf('%d of %d HRTFs equalized...\n',kk,size(samplingGrid,1))
    end
end

% Transform equalized HRTFs to SH-domain
fprintf('Performing spherical Fourier transform with N = %d. This may take some time...\n',N);
eqHRTFdataset = supdeq_hrtf2sfd(HRTF_equalized_L,HRTF_equalized_R,N,samplingGrid,fs,'ak',tikhEps); %Use AKtools...
eqHRTFdataset.HRTF_L = HRTF_equalized_L;
eqHRTFdataset.HRTF_R = HRTF_equalized_R;
eqHRTFdataset.FFToversize = sparseHRTFdataset.FFToversize;
eqHRTFdataset.eqEarDistance = eqDataset.earDistance;
eqHRTFdataset.eqWaveType = eqDataset.waveType;
if eqDataset.waveType == 1
    eqHRTFdataset.sourceDistance = eqDataset.sourceDistance;
end
if isfield(eqDataset,'limited')
    if eqDataset.limited
        eqHRTFdataset.limited = true;
        eqHRTFdataset.limitInfo = eqDataset.limitInfo;
    end
end
fprintf('Done with equalization...\n');

end

