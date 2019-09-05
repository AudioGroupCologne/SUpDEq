%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% script get_sparse_HRTF_set_Gauss
%
% Just a small script to show how to get a sparse HRTF dataset based on the
% reference HRTF dataset 'HRIRs_sfd_Nx' or 'HRIRs_NF150_sfd_Nx'.
% Field sourceDistance is used if available. Important for DVF (distance
% variation function) processing with 'HRIRs_NFxxx' as reference.
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

%%
clear all;
%Load reference dataset
referenceHRTFdataset = importdata('HRIRs_sfd_N35.mat');

%Simply do subsampling of the referenceHRTFdataset
%Define spatial order
nSparse = 1;
%Get sparse sampling grid
[sparseGrid, ~, NmaxSparse] = supdeq_gauss(nSparse);

%% Create HRTF dataset as mat file (for internal use)
%Get HRTFs from referenceHRTFdataset
[sparseHRTFdataset.HRTF_L, sparseHRTFdataset.HRTF_R] = supdeq_getArbHRTF(referenceHRTFdataset,sparseGrid);
%Fill struct with additional info
sparseHRTFdataset.f = referenceHRTFdataset.f;
sparseHRTFdataset.Nmax = NmaxSparse;
sparseHRTFdataset.FFToversize = referenceHRTFdataset.FFToversize;
sparseHRTFdataset.samplingGrid = sparseGrid;
if isfield(referenceHRTFdataset,'sourceDistance')
    sparseHRTFdataset.sourceDistance = referenceHRTFdataset.sourceDistance;
end

% Save as mat file
fileName = ['sparseHRTFdataset_G',num2str(nSparse)];
save(fileName,'sparseHRTFdataset')

%% Create HRIR dataset in SOFA format (for standardized input format)
%Get HRIRs from referenceHRTFdataset
[sparseHRIRdataset.HRIR_L, sparseHRIRdataset.HRIR_R] = supdeq_getArbHRIR(referenceHRTFdataset,sparseGrid);
% Write SOFA object
if isfield(referenceHRTFdataset,'sourceDistance')
    %Get source distance if is field
    sparseHRIRdataset_SOFA = supdeq_writeSOFAobj(sparseHRIRdataset.HRIR_L',sparseHRIRdataset.HRIR_R',sparseGrid(:,1:2),[],[],referenceHRTFdataset.sourceDistance);
else
    % %Use defaults: fs = 48000, earDistance = 0.165m, sourceDistance = 3.0m
    sparseHRIRdataset_SOFA = supdeq_writeSOFAobj(sparseHRIRdataset.HRIR_L',sparseHRIRdataset.HRIR_R',sparseGrid(:,1:2));
end

% Save as SOFA file
fileName = ['sparseHRIRdataset_G',num2str(nSparse),'.sofa'];
SOFAsave(fileName,sparseHRIRdataset_SOFA)
