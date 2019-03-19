%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% script get_HRIRs_sfd_N35
%
% Just a small script to show how to get the reference dataset
% 'HRIRs_sfd_N35' 
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

%%
% Load HRIR_L2702 SOFA dataset
hrirs = SOFAload('HRIR_L2702.sofa');

% Load respective spatial sampling grid. We use a grid with sampling
% weights in this case in order to use SOFiA toolbox for the transform. The
% sampling grid (and the corresponding weights) is derived from the 
% original miro file "HRIR_L2702.m" and transformed to degree.
samplingGrid = importdata('samplingGrid_HRIR_L2702.mat');

% Transform with order N = 35 and FFToversize = 4
HRIRs_sfd_N35 = supdeq_sofa2sfd(hrirs,35,samplingGrid,4);

% Save as mat file
save HRIRs_sfd_N35 HRIRs_sfd_N35;