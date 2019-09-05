%% SUpDEq - Spatial Upsampling by Directional Equalization
%  
% script supdeq_demo_ALFE
%
% This script shows how the ALFE (Adaptive Low Frequency Extention)
% algorithm can be applied in combination with SUpDEq.
%
% Reference:
% C. Pörschmann and J. M. Arend, "Obtaining Dense HRTF Sets from Sparse 
% Measurements in Reverberant Environments," in Proceedings of the 
% AES International Conference on Immersive and Interactive Audio, 
% York, UK, 2019, pp. 1?10.
%
% (C) 2019 by JMA, Johannes M. Arend
%             CP, Christoph Pörschmann 
%             Technische Hochschule Köln
%             University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing
%
%% (1) - Load reference HRTF dataset (stored as SH-coefficients in spherical harmonics domain)
%See publication for more info about the measurement setup of the reference HRTF dataset.
%The reference dataset is only needed for comparison with the final de-equalized dataset...
referenceHRTFdataset = importdata('RAR_HRIR_1m_Genelec_2702SFD.mat');

%% (2) - Load sparse HRTF dataset (stored as FFT-bins in Fourier domain)
%See publication for more info about the measurement setup of the sparse HRTF dataset.
load('ZW8_3_HRIR_1m_Genelec_sp_L86.mat');

%% (3) - Get equalization dataset (SH-coefficients)
%The eqDataset describes the sound pressure distribution on a sphere 
%Use defaults: N = 35, earDistance = 0.165m, NFFT = 512, fs = 48000;
NFFT=(length(sparseHRTFdataset.HRTF_L(1,:))-1)*2;
eqDataset = supdeq_getEqDataset(35,0.165,NFFT,48000);

%% (4) - Perform equalization
%Here, the sparse HRTF dataset is equalized with the eqDataset. The
%equalized HRTF are transformed to the spherical harmonics domain again with
%the maximal order N which is possible with the sparse sampling grid.
%N and the sparse sampling grid are part of the sparseHRTFdataset struct
sparseSamplingGrid = sparseHRTFdataset.samplingGrid;
Nsparse = sparseHRTFdataset.Nmax;

eqHRTFdataset = supdeq_eq(sparseHRTFdataset,eqDataset,Nsparse,sparseSamplingGrid);

%% (4a) - Perform ALFE - If desired
%In this section the lower frequency components below e certain frequency
%xover are replaced by a frequecy-shaped dirac. The ALFE algorithm is done
%on the equalized dataset in order to benefit from the constant frequency
%response towards low frequencies. 
xover = [500 750];
fs = 48000;
[eqHRTFdataset_ALFE] = supdeq_alfe(eqHRTFdataset,sparseHRTFdataset.samplingGrid,xover,fs);
eqHRTFdataset=eqHRTFdataset_ALFE;  

%% (5) - Perform de-equalization 
%Here, the sparse equalized HRTF dataset is de-equalized with the
%deqDataset. This is done on a dense spatial sampling grid. The results is a
%dense HRTF/HRIR dataset. In this case, deqDataset is the same as the
%eqDataset...

%First, define dense spatial sampling grid. Here, we use the lebedev grid
%with 2702 points again (same as the reference HRIR dataset).
%The highest stable grid order here is N = 44. However, we use N = 35 for the
%spherical Fourier transform of the de-equalized HRTFs.
denseSamplingGrid = supdeq_lebedev(2702);
Ndense = 35;

%Perform de-equalization. Apply head and tail window (8 and 32 samples
%respectively) to de-equalized HRIRs/HRTFs.
[denseHRTFdataset, denseHRIRdataset, denseHRTFdataset_sh] = supdeq_deq(eqHRTFdataset, eqDataset, Ndense, denseSamplingGrid,[8,32]);

%% (6) - Optional: Save as SOFA object
%Use defaults: fs = 48000, earDistance = 0.165m, sourceDistance = 3.0m
%denseHRIRdataset_SOFA = supdeq_writeSOFAobj(denseHRIRdataset.HRIR_L,denseHRIRdataset.HRIR_R,denseSamplingGrid);

%% (6a) - Optional: Save SH coefficients
save ('denseHRTFdataset_JBL_L2702.mat', 'denseHRTFdataset_sh')

%% (7) - Optional: Plot HRIRs
%Get HRIRs from reference dataset and de-equalized dense dataset. In this
%example, we chose a lateral source, because most differences between the
%reference and the de-equalized HRIRs occure at the contralateral ear.
azPlot = 180;
elPlot = 90;

[hrir_ref(:,1),hrir_ref(:,2)] = supdeq_getArbHRIR(referenceHRTFdataset,[azPlot,elPlot],'DEG',2,'ak');
[hrir_deq(:,1),hrir_deq(:,2)] = supdeq_getArbHRIR(denseHRTFdataset_sh,[azPlot,elPlot],'DEG',2,'ak');

supdeq_plotIR(hrir_ref(:,2),hrir_deq(:,2),[],[],1);

%% (8) - Optional: Listen to the results
%Here, you can listen to a short drums sequence, spatialized with a
%reference HRIR and a de-equalized HRIR

%Load drums test signal
[testSignal,fsTestSignal] = audioread('drums.wav');

%Define az and el for playback
azPlayback = 270;
elPlayback = 90;

%Convolve with reference HRIR and play back
supdeq_listen(referenceHRTFdataset,testSignal,[azPlayback,elPlayback]);
%Wait for length of test-signal
pause(length(testSignal)/fsTestSignal+0.1);
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
