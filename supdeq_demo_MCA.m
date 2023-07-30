%% SUpDEq - Spatial Upsampling by Directional Equalization
%  
% script supdeq_demo_MCA
%
% This script presents magnitude-corrected and time-aligned interpolation,
% hereafter abbreviated as MCA interpolation. The method combines a
% perceptually motivated post-interpolation magnitude correction with
% time-aligned interpolation. This demo script uses SUpDEq processing for
% time alignment and spherical harmonics (SH) interpolation for spatial
% upsampling of a sparse HRTF set to a dense HRTF set. The script compares
% the results of MCA interpolation and conventional time-aligned
% interpolation. The proposed magnitude correction as well as various 
% time-alignment and interpolation approaches are implemented in the
% function supdeq_interpHRTF, used throughout this demo script.
%
% Reference:
% J. M. Arend, C. Pörschmann, S. Weinzierl, F. Brinkmann, 
% "Magnitude-Corrected and Time-Aligned Interpolation of 
% Head-Related Transfer Functions," (Manuscript submitted for publication).
%
% (C) 2022/2023 by JMA, Johannes M. Arend
%               TU Berlin
%               Audio Communication Group

%% (1) - Define sparse and dense grid

%Sparse Lebedev Grid
%Azimuth, Colatitude, Weight
Ns = 3;
sgS = supdeq_lebedev([],Ns);

%Dense Lebedev Grid (Target grid for upsampling)
%Azimuth, Colatitude, Weight
Nd = 44;
sgD = supdeq_lebedev([],Nd);

%% (2) - Get sparse HRTF set

%Get sparse KU100 HRTF set based on "HRIR_L2702.sofa" dataset in SH domain

sparseHRTF = supdeq_getSparseDataset(sgS,Ns,44,'ku100');
fs = sparseHRTF.fs;

%% (3) - Perform spatial upsampling to dense grid

headRadius = 0.0875; %Also default value in function

%Interpolation / upsampling with SH interpolation but without any
%pre/post-processing ('none'), called SH here
interpHRTF_sh = supdeq_interpHRTF(sparseHRTF,sgD,'None','SH',nan,headRadius);

%Interpolation / upsampling with SUpDEq time alignment and SH
%interpolation, called 'conventional'
interpHRTF_con = supdeq_interpHRTF(sparseHRTF,sgD,'SUpDEq','SH',nan,headRadius);

%Interpolation / upsampling with MCA interpolation, i.e., in this example
%SUpDEq time alignment and SH interpolation plus post-interpolation
%magnitude correction, as described in the paper
%Maximum boost of magnitude correction set to inf (unlimited gain), 
%resulting in no soft-limiting (as presented in the paper)
%The other MCA settings are by default as described in the paper, i.e.:
% (a) Magnitude-correction filters are set to 0dB below spatial aliasing
% frequency fA
% (b) Soft-limiting knee set to 0dB (no knee)
% (c) Magnitude-correction filters are designed as minimum phase filters
interpHRTF_mca = supdeq_interpHRTF(sparseHRTF,sgD,'SUpDEq','SH',inf,headRadius);

%The resulting datasets can be saved as SOFA files using the function
%supdeq_writeSOFAobj

%% (4) - Optional: Get reference HRTF set

%Get reference dataset (using the "supdeq_getSparseDataset" function for
%convenience here)
refHRTF = supdeq_getSparseDataset(sgD,Nd,44,'ku100');

%Also get HRIRs
refHRTF.HRIR_L = ifft(AKsingle2bothSidedSpectrum(refHRTF.HRTF_L.'));
refHRTF.HRIR_L = refHRTF.HRIR_L(1:end/refHRTF.FFToversize,:);
refHRTF.HRIR_R = ifft(AKsingle2bothSidedSpectrum(refHRTF.HRTF_R.'));
refHRTF.HRIR_R = refHRTF.HRIR_R(1:end/refHRTF.FFToversize,:);

%% (5) - Optional: Plot HRIRs

%Plot frontal left-ear HRIR (Az = 0, Col = 90) of conventional and MCA
%interpolated dataset. Differences are small.
idFL = 16; %ID in sgD
supdeq_plotIR(interpHRTF_con.HRIR_L(:,idFL),interpHRTF_mca.HRIR_L(:,idFL),[],[],8);

%Plot contralateral left-ear HRIR (Az = 270, Col = 90) of SH-only and MCA
%interpolated dataset. Differences are strong
idCL = 2042; %ID in sgD
supdeq_plotIR(interpHRTF_sh.HRIR_L(:,idCL),interpHRTF_mca.HRIR_L(:,idCL),[],[],8);

%Plot contralateral left-ear HRIR (Az = 270, Col = 90) of conventional and MCA
%interpolated dataset. Differences are still significant.
supdeq_plotIR(interpHRTF_con.HRIR_L(:,idCL),interpHRTF_mca.HRIR_L(:,idCL),[],[],8);

%Plot contralateral left-ear HRIR (Az = 270, Col = 90) of SH-only 
%interpolated dataset and reference.
supdeq_plotIR(interpHRTF_sh.HRIR_L(:,idCL),refHRTF.HRIR_L(:,idCL),[],[],8);

%Plot contralateral left-ear HRIR (Az = 270, Col = 90) of conventional 
%interpolated dataset and reference.
supdeq_plotIR(interpHRTF_con.HRIR_L(:,idCL),refHRTF.HRIR_L(:,idCL),[],[],8);

%Plot contralateral left-ear HRIR (Az = 270, Col = 90) of MCA 
%interpolated dataset and reference.
%MCA interpolated HRTF is much closer to the reference than HRTF from
%conventional interpolation. The "bump" above 10 kHz is corrected through
%the magnitude correction. Interpolation errors (in this case spatial
%aliasing errors at the contralateral ear) are reduced
supdeq_plotIR(interpHRTF_mca.HRIR_L(:,idCL),refHRTF.HRIR_L(:,idCL),[],[],8);

%% (6) - Optional: Calculate log-spectral difference 

%Calculate left-ear log-spectral differences in dB.
%Log-spectral difference is often used as a measure for spectral distance
%between a reference and an interpolated HRTF set. Not the same as 
%the magnitude error in auditory filters presented in the paper!

lsd_sh =  supdeq_calcLSD_HRIR(interpHRTF_sh.HRIR_L,refHRTF.HRIR_L,fs,16);
lsd_con = supdeq_calcLSD_HRIR(interpHRTF_con.HRIR_L,refHRTF.HRIR_L,fs,16);
lsd_mca = supdeq_calcLSD_HRIR(interpHRTF_mca.HRIR_L,refHRTF.HRIR_L,fs,16);

%Quick plot of log-spectral difference averaged over all directions of the
%2702-point Lebedev grid. As shown in the paper, MCA provides the most
%significant benefit in the critical contralateral region. Thus, when 
%averaging over all positions, the benefit of MCA seems smaller. The
%analysis in the paper as well as the provided audio examples reveal
%in much more detail the considerable improvements that the proposed 
%magnitude correction provides.
AKf(18,9);
semilogx(lsd_sh.f,lsd_sh.lsd_freq,'LineWidth',1.5);
hold on;
semilogx(lsd_con.f,lsd_con.lsd_freq,'LineWidth',1.5);
semilogx(lsd_mca.f,lsd_mca.lsd_freq,'LineWidth',1.5);
xlim([500 20000])
legend('SH W/O MC','SH SUpDEq W/O MC','SH SUpDEq W/ MC','Location','NorthWest');
xlabel('Frequency in Hz');
ylabel('Log-Spectral Difference in dB');
grid on;

%% (7) - Optional: Calculate magnitude error in auditory filters (similar to paper)

[erb_sh,fc_erb] = AKerbError(interpHRTF_sh.HRIR_L,refHRTF.HRIR_L, [50 fs/2], fs);
erb_con = AKerbError(interpHRTF_con.HRIR_L,refHRTF.HRIR_L, [50 fs/2], fs);
erb_mca = AKerbError(interpHRTF_mca.HRIR_L,refHRTF.HRIR_L, [50 fs/2], fs);

%Quick plot of absolute magnitude error in auditory filters averaged over all
%directions of the 2702-point Lebedev grid. 
AKf(18,9);
semilogx(fc_erb,mean(abs(erb_sh),2),'LineWidth',1.5);
hold on;
semilogx(fc_erb,mean(abs(erb_con),2),'LineWidth',1.5);
semilogx(fc_erb,mean(abs(erb_mca),2),'LineWidth',1.5);
xlim([500 20000])
legend('SH W/O MC','SH SUpDEq W/O MC','SH SUpDEq W/ MC','Location','NorthWest');
xlabel('Frequency in Hz');
ylabel('\DeltaG(f_c) in dB');
grid on;





