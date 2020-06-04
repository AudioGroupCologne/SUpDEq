%% SUpDEq - Spatial Upsampling by Directional Equalization
%  
% script supdeq_demo_DVF
%
% This script presents how SUpDEq can be applied as a distance variation
% function (DVF). A (sparse) far-field HRTF set is transformed to (dense)
% near-field HRTF set. 
%
% Reference:
% J. M. Arend and C. Pörschmann, "Synthesis of Near-Field HRTFs by Directional 
% Equalization of Far-Field Datasets," in Proceedings of the 45th DAGA, 2019, pp. 1454?1457.
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

%% (1) - Load reference HRTF dataset (stored as SH-coefficients in spherical harmonics domain)
%See small script "get_HRIRs_sfd_Nx.m" in materials folder to see how the
%dataset was obtained based on the "HRIR_L2702_NFxxx.sofa" HRIR dataset of a
%Neumann KU100 dummy head. The reference dataset is only needed for
%comparison with the final de-equalized dataset...

referenceHRTFdataset = importdata('HRIRs_NF050_sfd_N44.mat');
NeqDataset = 44; %According to referenceHRTFdataset

%% (2) - Load sparse HRTF dataset (stored as FFT-bins in Fourier domain)
%This is just an example. You could use any sparse HRTF dataset here! Also
%have a look at the script "get_sparse_HRTF_set.m" in materials folder to
%see how the sparse datasets are obtained.

%Here, we use a sparse HRTF daset with only 38 sample points (Lebedev grid)
%The sampling grid as well as the highste possible transform oder N are
%required for the equalization! It is saved in the struct here, but it
%could also be passed separately to the function "subdeq_eq"
load('sparseHRTFdataset_DVF_L38.mat');

%Optionally, a sparse HRIR dataset in SOFA format can be loaded. With the
%function supdeq_sofa2hrtf, the SOFA file will be transformed to
%frequency domain and saved in a struct as used by the SUpDEq toolbox
%(compare to the structs 'sparseHRTFdataset'). FFToversize and Nmax have to
%be defined manually now, because SOFA does not support such meta data.
%Furthermore, the sampling grid in SOFA does not support sampling weights,
%which can be useful for spherical Fourier / SH transform. Thus, optionally
%a sampling grid with weigths can be passed, which will than be stored in
%the struct. Otherwise, the sampling grid stored in the SOFA file will be 
%written to the HRTF struct.
SOFAinput = false;
if SOFAinput
    %Load sparse HRIR dataset in SOFA format
    sparseHRIRdataset_SOFA = SOFAload('sparseHRIRdataset_DVF_L38.sofa');
    %Transform to sparseHRTFdataset struct with pre-defined samplingGrid 
    %(Lebedev grid with 38 nodes here), Nmax = 4, and FFToversize = 4.
    sparseHRTFdataset = supdeq_sofa2hrtf(sparseHRIRdataset_SOFA,4,supdeq_lebedev(38),4);
end

%% (3) - Get equalization dataset (SH-coefficients)
%The eqDataset describes the sound pressure distribution on a sphere 
%Use N = 44 and defaults: earDistance = 0.165m, NFFT = 512, fs = 48000;
eqDataset = supdeq_getEqDataset(NeqDataset);

%% (4) - Perform equalization
%Here, the sparse HRTF dataset is equalized with the eqDataset. The
%equalized HRTF are transformed to the SH-domain again with the maximal 
%order N which is possible with the sparse sampling grid.
%N and the sparse sampling grid are part of the sparseHRTFdataset struct
sparseSamplingGrid = sparseHRTFdataset.samplingGrid;
Nsparse = sparseHRTFdataset.Nmax;

eqHRTFdataset = supdeq_eq(sparseHRTFdataset,eqDataset,Nsparse,sparseSamplingGrid);

%% (5) - Perform de-equalization 
%Here, the sparse equalized HRTF dataset is de-equalized with the
%deqDataset. This is done on a dense spatial sampling grid. The results is a
%dense HRTF/HRIR dataset. To adapt the distance (apply the DVF), the
%deqDataset is generated based on a spherical wave at the target distance

%First, define dense spatial sampling grid. Here, we use the lebedev grid
%with 2702 points again (same as the reference HRIR dataset).
%The highest stable grid order here is N = 44.
denseSamplingGrid = supdeq_lebedev(2702);
Ndense = 44;

%Now, generate deqDataset for targetDistance
targetDistance = 0.50; %Target distance for the near-field HRTFs in m
%Use defaults: earDistance = 0.165m, NFFT = 512, fs = 48000
deqDataset = supdeq_getEqDataset(NeqDataset,[],[],[],1,targetDistance);

%Perform de-equalization. Apply head and tail window (8 and 32 samples
%respectively) to de-equalized HRIRs/HRTFs. The directions of the
%denseSamplingGrid are re-calculated automatically in supdeq_deq
%to consider the acoustic parallax effect which occurs in the near field
[denseHRTFdataset, denseHRIRdataset, denseHRTFdataset_sh] = supdeq_deq(eqHRTFdataset, deqDataset, Ndense, denseSamplingGrid,[8,32],0.99);

%% (6) - Optional: Save as SOFA object
%Use sourceDistance = targetDistance and defaults: fs = 48000, earDistance = 0.165m
denseHRIRdataset_SOFA = supdeq_writeSOFAobj(denseHRIRdataset.HRIR_L,denseHRIRdataset.HRIR_R,denseSamplingGrid,[],[],targetDistance);

%% (7) - Optional: Plot HRIRs
%Get HRIRs from reference dataset and de-equalized dense dataset. In this
%example, we chose a lateral source, because most differences between the
%reference and the de-equalized HRIRs occure at the contralateral ear.
azPlot = 0;
elPlot = 90;
[hrir_ref(:,1),hrir_ref(:,2)] = supdeq_getArbHRIR(referenceHRTFdataset,[azPlot,elPlot]);
[hrir_deq(:,1),hrir_deq(:,2)] = supdeq_getArbHRIR(denseHRTFdataset_sh,[azPlot,elPlot]);

%Plot left and right channel of reference and deq set in two plots with 
%32 x FFT oversampling
supdeq_plotIR(hrir_ref(:,1),hrir_deq(:,1),[],[],32);
supdeq_plotIR(hrir_ref(:,2),hrir_deq(:,2),[],[],32);

%% (8) - Optional: Listen to the results
%Here, you can listen to a short drums sequence, spatialized with a
%reference HRIR and a de-equalized HRIR

%Load drums test signal
[testSignal,fsTestSignal] = audioread('drums.wav');

%Define az and el for playback
azPlayback = 0;
elPlayback = 90;

%Convolve with reference HRIR and play back
supdeq_listen(referenceHRTFdataset,testSignal,[azPlayback,elPlayback]);
%Wait for length of test-signal
pause(length(testSignal)/fsTestSignal+0.05);
%Convolve with de-equalized HRIR and play back
supdeq_listen(denseHRTFdataset_sh,testSignal,[azPlayback,elPlayback]);

%% (9) - Optional: Technical Evaluation

%Perform test of sagittal plane localization with default spatial sampling
%grid
spLocalization = supdeq_testSPlocalization(denseHRTFdataset_sh, referenceHRTFdataset);

%Perform test of horizontal plane localizatoin with default spatial
%sampling grid
hpLocalization = supdeq_testHPlocalization(denseHRTFdataset_sh, referenceHRTFdataset);

%Calculate spectral differences with default spatial sampling grid and 
%make simple plot
specDifferences = supdeq_calcSpectralDifference(denseHRTFdataset_sh, referenceHRTFdataset,[],[],false,true);

%Calculate spectral differences in auditory bands with default spatial 
%sampling grid and make simple plot
specDifferences_auditoryBands = supdeq_calcSpectralDifference(denseHRTFdataset_sh, referenceHRTFdataset,[],[],true,true);

%Calculate energy per order N and make simple plot
energyPerOrder = supdeq_calcEnergyPerOrder(denseHRTFdataset_sh,true);
