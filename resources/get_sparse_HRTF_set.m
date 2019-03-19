%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% script get_sparse_HRTF_set
%
% Just a small script to show how to get a sparse HRTF dataset based on the
% reference HRTF dataset 'HRIRs_sfd_N35' 
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
%Define number of nodes of sparse grid. Options for lebedev grid:
%6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,
%350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,
nSparse=86;
%Get sparse sampling grid
[sparseGrid, ~, NmaxSparse] = supdeq_lebedev(nSparse);

%% Create HRTF dataset as mat file (for internal use)
%Get HRTFs from referenceHRTFdataset
[sparseHRTFdataset.HRTF_L, sparseHRTFdataset.HRTF_R] = supdeq_getArbHRTF(referenceHRTFdataset,sparseGrid);
%Fill struct with additional info
sparseHRTFdataset.f = referenceHRTFdataset.f;
sparseHRTFdataset.Nmax = NmaxSparse;
sparseHRTFdataset.FFToversize = referenceHRTFdataset.FFToversize;
sparseHRTFdataset.samplingGrid = sparseGrid;

% Save as mat file
fileName = ['sparseHRTFdataset_L',num2str(nSparse)];
save(fileName,'sparseHRTFdataset')

%% Create HRIR dataset in SOFA format (for standardized input format)
%Get HRIRs from referenceHRTFdataset
[sparseHRIRdataset.HRIR_L, sparseHRIRdataset.HRIR_R] = supdeq_getArbHRIR(referenceHRTFdataset,sparseGrid);
% Write SOFA object
% %Use defaults: fs = 48000, earDistance = 0.165m, sourceDistance = 3.0m
sparseHRIRdataset_SOFA = supdeq_writeSOFAobj(sparseHRIRdataset.HRIR_L',sparseHRIRdataset.HRIR_R',sparseGrid(:,1:2));

% Save as SOFA file
fileName = ['sparseHRIRdataset_L',num2str(nSparse),'.sofa'];
SOFAsave(fileName,sparseHRIRdataset_SOFA)
