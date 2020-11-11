%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% script supdeq_compareTimeAlignmentMethods
% 
% In this script, various preprocessing / time-alignment methods [1,2,3,5] are implemented
% and compared to the SUpDEq method [4]. It is the basis for the journal submission [6].
% In brief, the sparse HRTF set gets time-aligned with the respective method, interpolated
% to a dense HRTF set in SH domain, and then HRTFs get reconstructed with the respective
% inverse alignment. The methods are described and compared in more detail in [6].
%
% Dependencies: SUpDEq toolbox with SOFiA and AKtools
%               https://github.com/AudioGroupCologne/SUpDEq
%
% References:
% Onset-Based Time-Alignment (OBTA)
% [1] M. J. Evans, J. A. S. Angus, and A. I. Tew, "Analyzing head-related transfer function measurements 
% using surface spherical harmonics," 
% J. Acoust. Soc. Am., vol. 104, no. 4, pp. 2400-2411, 1998.
% [2] F. Brinkmann and S. Weinzierl, "Comparison of head-related transfer functions pre-processing 
% techniques for spherical harmonics decomposition," 
% in Proceedings of the AES Conference on Audio for Virtual and Augmented Reality, 2018, pp. 1-10.
%
% Frequency-Dependenten Time-Alignment (FDTA)
% [3] M. Zaunschirm, C. Schoerkhuber, and R. Hoeldrich, 
% "Binaural rendering of Ambisonic signals by HRIR time alignment and a diffuseness constraint," 
% J. Acoust. Soc. Am., vol. 143, no. 6, pp. 3616-3627, 2018.
%
% Spatial Upsampling by Directional Equalization (SUpDEq)
% [4] C. Pörschmann, J. M. Arend, and F. Brinkmann, 
% "Directional Equalization of Sparse Head-Related Transfer Function Sets for Spatial Upsampling," 
% IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 27, no. 6, pp. 1060-1071, 2019.
%
% Phase-Correction (PC)
% [5] Z. Ben-Hur, D. Lou Alon, R. Mehra, and B. Rafaely, 
% "Efficient Representation and Sparse Sampling of Head-Related Transfer Functions 
% Using Phase-Correction Based on Ear Alignment," 
% IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 27, no. 12, pp. 2249-2262, 2019.
%
% Submission
% [6] J.M. Arend, F. Brinkmann, and C. Pörschmann, 
% "Assessing Spherical Harmonics Interpolation of Time-Aligned Head-Related
% Transfer Functions"
% Manuscript submitted for publication, 2020.
%
% (C) 2020 by JMA, Johannes M. Arend
%             Technische Hochschule Köln
%             University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing
%
%%
clear all;
close all;

%Chose SH order of sparse HRTF dataset (Lebedev sampling scheme)
sparseOrder = 4;
%Get sparse HRTF set
[sparseHRTFdataset,sparseHRIRdataset] = supdeq_getSparseDataset(supdeq_lebedev([],sparseOrder),sparseOrder,44);
%Get fs
fs = sparseHRTFdataset.f(end)*2;

%Define dense sampling grid for upsampling
denseSamplingGrid = supdeq_lebedev(2702);
%Load reference dataset from HD
referenceHRTFdataset = importdata('HRIRs_sfd_N44.mat');
%Get Ndense, which in this case is the N = 44 (order of
%supdeq_lebedev(2702)
Ndense = referenceHRTFdataset.N;

%Calculate optimal radius for KU100 according to Algazi et al. (2001)
rKU100 = supdeq_optRadius(0.155,0.25,0.20,'Algazi'); 
%Define c
c = 343;

%% Perform Onset-Based Time-Alignment (OBTA)
%Implementation according to [2]

%Estimate TOA on low-passed and 10 times upsampled HRIRs according to [2]
toa_L = AKonsetDetect(sparseHRIRdataset.HRIR_L, 10, -20, 'rel', [3e3 fs]);
toa_R = AKonsetDetect(sparseHRIRdataset.HRIR_R, 10, -20, 'rel', [3e3 fs]);

%Transform TOAs to SH domain
toa_L_nm = AKsht(toa_L, false, sparseHRIRdataset.samplingGrid, sparseHRIRdataset.Nmax, 'complex', fs, false, 'real');
toa_R_nm = AKsht(toa_R, false, sparseHRIRdataset.samplingGrid, sparseHRIRdataset.Nmax, 'complex', fs, false, 'real');

%Time-Align HRIRs (remove TOAs) with fractional delay
HRIR_L_OBTA = AKfractionalDelayCyclic(sparseHRIRdataset.HRIR_L, -toa_L);
HRIR_R_OBTA = AKfractionalDelayCyclic(sparseHRIRdataset.HRIR_R, -toa_R);

%Transform to SH domain with Nmax according to sparse dataset
HRTF_L_OBTA = AKboth2singleSidedSpectrum (fft(HRIR_L_OBTA));
HRTF_R_OBTA = AKboth2singleSidedSpectrum (fft(HRIR_R_OBTA));
obtaHRTFdataset = supdeq_hrtf2sfd(HRTF_L_OBTA.',HRTF_R_OBTA.',sparseHRTFdataset.Nmax,sparseHRTFdataset.samplingGrid,fs,'ak');

%Upsample to dense grid by inverse SH transform
[HRIR_L_Dense_OBTA, HRIR_R_Dense_OBTA] = supdeq_getArbHRIR(obtaHRTFdataset,denseSamplingGrid,'DEG',2,'ak');

%Get interpolated TOAs by inverse SH transform
toa_L_Dense = AKisht(toa_L_nm, false, denseSamplingGrid(:,1:2), 'complex', false, false, 'real');
toa_R_Dense = AKisht(toa_R_nm, false, denseSamplingGrid(:,1:2), 'complex', false, false, 'real');

%Re-Insert TOAs with fractional delays
HRIR_L_Dense_OBTA = AKfractionalDelayCyclic(HRIR_L_Dense_OBTA, toa_L_Dense);
HRIR_R_Dense_OBTA = AKfractionalDelayCyclic(HRIR_R_Dense_OBTA, toa_R_Dense);

%Transform dense HRIR dataset to SH domain (apply zero padding according to
% reference) at high order Ndense
HRTF_L_Dense_OBTA = AKboth2singleSidedSpectrum (fft(HRIR_L_Dense_OBTA,size(HRIR_L_Dense_OBTA,1)*referenceHRTFdataset.FFToversize));
HRTF_R_Dense_OBTA = AKboth2singleSidedSpectrum (fft(HRIR_R_Dense_OBTA,size(HRIR_R_Dense_OBTA,1)*referenceHRTFdataset.FFToversize));
denseHRTFdataset_obta = supdeq_hrtf2sfd(HRTF_L_Dense_OBTA.',HRTF_R_Dense_OBTA.',Ndense,denseSamplingGrid,fs,'ak');
denseHRTFdataset_obta.FFToversize = referenceHRTFdataset.FFToversize;

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%% Perform Frequency-Dependente Time-Alignment (FDTA)

%Transform sampling grid as required with azimuth between -pi to pi and
%elevation -pi/2 to pi/2
nPoints = size(sparseHRTFdataset.samplingGrid,1);
az = sparseHRTFdataset.samplingGrid(:,1);
az(az>180) = (360-az(az>180))*-1;
el = 90-sparseHRTFdataset.samplingGrid(:,2);
azRad = az * pi / 180;
elRad = el * pi / 180;

% Calculate the relative TOAs for each grid point according to Eq. 16 in [3]
tl = cos(elRad).*sin(azRad)*rKU100*(c^-1); 
tr = -tl;

%Get normalized cut on frequency at f = 1500Hz
fVec =  sparseHRTFdataset.f;
fCut = 1500;
[~,fCutBin] = min(abs(fVec-fCut));

%Design allpass filter for each grid point according to Eq.26 in [3]
wc = 2*pi*fVec(fCutBin);
apl = zeros(nPoints,length(fVec));
apl(:,1:fCutBin-1) = 1;
apr = zeros(nPoints,length(fVec));
apr(:,1:fCutBin-1) = 1;
for id = 1:nPoints
    for kk = fCutBin:size(apl,2)
        apl(id,kk) = exp(-1j*(2*pi*fVec(kk)-wc)*tl(id));
        apr(id,kk) = exp(-1j*(2*pi*fVec(kk)-wc)*tr(id));
    end
end

%Apply allpass to HRTFs
HRTF_L_FDTA = sparseHRTFdataset.HRTF_L .* apl;
HRTF_R_TAZ = sparseHRTFdataset.HRTF_R .* apr;

%Transform to SH domain with Nmax according to sparse dataset
fdtaHRTFdataset = supdeq_hrtf2sfd(HRTF_L_FDTA,HRTF_R_TAZ,sparseHRTFdataset.Nmax,sparseHRTFdataset.samplingGrid,fs,'ak');

%Upsample to dense grid by inverse SH transform
[HRTF_L_Dense_FDTA, HRTF_R_Dense_FDTA] = supdeq_getArbHRTF(fdtaHRTFdataset,denseSamplingGrid,'DEG',2,'ak');

%Add relative TOAs for dense sampling grid to interpolated HRTFs
nPointsDense = size(denseSamplingGrid,1);
azDense = denseSamplingGrid(:,1);
azDense(azDense>180) = (360-azDense(azDense>180))*-1;
elDense = 90-denseSamplingGrid(:,2);
azDenseRad = azDense * pi / 180;
elDenseRad = elDense * pi / 180;

%Design allpass filter for each grid point of dense grid according to Eq.26 in [3]
tlDense = cos(elDenseRad).*sin(azDenseRad)*rKU100*(c^-1);
trDense = -tlDense;

aplDense = zeros(nPointsDense,length(fVec));
aplDense(:,1:fCutBin-1) = 1;
aprDense = zeros(nPointsDense,length(fVec));
aprDense(:,1:fCutBin-1) = 1;
for id = 1:nPointsDense
    for kk = fCutBin:size(apl,2)
        %No minus before 1j to get inverse phase term. Could also be achieved by division instead of multiplication
        aplDense(id,kk) = exp(1j*(2*pi*fVec(kk)-wc)*tlDense(id)); 
        aprDense(id,kk) = exp(1j*(2*pi*fVec(kk)-wc)*trDense(id));
    end
end

%Add relative TOAs of dense sampling grid to interpolated HRTFs
HRTF_L_Dense_FDTA = HRTF_L_Dense_FDTA .* aplDense;
HRTF_R_Dense_FDTA = HRTF_R_Dense_FDTA .* aprDense;

%Transform dense HRTF dataset to SH domain 
denseHRTFdataset_fdta = supdeq_hrtf2sfd(HRTF_L_Dense_FDTA,HRTF_R_Dense_FDTA,Ndense,denseSamplingGrid,fs,'ak');
denseHRTFdataset_fdta.FFToversize = referenceHRTFdataset.FFToversize;

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_fdta denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%% Perform Spatial Upsampling by Directional Equalization (SUpDEq)

%Get eqDataset
eqDataset = supdeq_getEqDataset(Ndense,2*rKU100,length(referenceHRTFdataset.f)*2-2,fs);

%Limit eqDataset (Small improvement explained in [6] leading to better
%results below the spatial aliasing frequency fA)
eqDataset = supdeq_limitEqDataset(eqDataset,sparseHRTFdataset.Nmax,eqDataset.radius);

%Equalization
[eqHRTFdataset, HRTF_equalized_L, HRTF_equalized_R] = supdeq_eq(sparseHRTFdataset,eqDataset,sparseHRTFdataset.Nmax,sparseHRTFdataset.samplingGrid);

%De-Equalization
[~,~,denseHRTFdataset_supdeq] = supdeq_deq(eqHRTFdataset, eqDataset, Ndense, denseSamplingGrid);

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_supdeq denseHRTFdataset_fdta denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%% Perform Phase-Correction (PC)

%Get wavenumber k
k  = 2*pi*sparseHRTFdataset.f/c; k = k.';

%Transform sampling grid to radiant
sg = sparseHRTFdataset.samplingGrid * pi / 180;

%Get phase correction term for left/right ear according to Eq. 13 in [5]
cosThetaL = cos(sg(:,2)')*cos(pi/2) + sin(sg(:,2)')*sin(pi/2) .* cos(sg(:,1)'-pi/2); %Left ear with -pi/2
phaseCorL = exp(-1j*rKU100 * k .* cosThetaL);
cosThetaR = cos(sg(:,2)')*cos(pi/2) + sin(sg(:,2)')*sin(pi/2) .* cos(sg(:,1)'+pi/2); %Right ear with +pi/2
phaseCorR = exp(-1j*rKU100 * k .* cosThetaR);

%Apply to HRTFs
HRTF_L_PC = sparseHRTFdataset.HRTF_L .* phaseCorL.'; %Just flip, no conjugate complex
HRTF_R_PC = sparseHRTFdataset.HRTF_R .* phaseCorR.';

%Transform to SH domain with Nmax according to sparse dataset
pcHRTFdataset = supdeq_hrtf2sfd(HRTF_L_PC,HRTF_R_PC,sparseHRTFdataset.Nmax,sparseHRTFdataset.samplingGrid,fs,'ak');

%Upsample to dense grid by inverse SH transform
[HRTF_L_Dense_PC, HRTF_R_Dense_PC] = supdeq_getArbHRTF(pcHRTFdataset,denseSamplingGrid,'DEG',2,'ak');

%Transform dense sampling grid to radiant
denseSamplingGridRad = denseSamplingGrid * pi / 180;

%Add relative TOAs of dense sampling grid to HRTFs (inverse phase
%correction)
cosThetaDenseL = cos(denseSamplingGridRad(:,2)')*cos(pi/2) + sin(denseSamplingGridRad(:,2)')*sin(pi/2) .* cos(denseSamplingGridRad(:,1)'+pi/2); %Inverse phase has +- pi at end of cosTheta equation
%No minus before 1j to get inverse phase term. Could also be achieved by division instead of multiplication
phaseCorDenseL = exp(1j*rKU100 * k .* cosThetaDenseL);
cosThetaDenseR = cos(denseSamplingGridRad(:,2)')*cos(pi/2) + sin(denseSamplingGridRad(:,2)')*sin(pi/2) .* cos(denseSamplingGridRad(:,1)'-pi/2);
%No minus before 1j to get inverse phase term. Could also be achieved by division instead of multiplication
phaseCorDenseR = exp(1j*rKU100 * k .* cosThetaDenseR);
HRTF_L_Dense_PC = HRTF_L_Dense_PC .* phaseCorDenseL.'; %Just flip, no conjugate complex
HRTF_R_Dense_PC = HRTF_R_Dense_PC .* phaseCorDenseR.';

%Transform dense HRTF dataset to SH domain
denseHRTFdataset_pc = supdeq_hrtf2sfd(HRTF_L_Dense_PC,HRTF_R_Dense_PC,Ndense,denseSamplingGrid,fs,'ak');
denseHRTFdataset_pc.FFToversize = referenceHRTFdataset.FFToversize;

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_pc denseHRTFdataset_supdeq denseHRTFdataset_fdta denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Some plots for physical evaluation. See more detailed results in [6] and
%in supplementary material of [6].
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Spectral difference over the entire sphere (2702 Lebedev grid)

supdeq_calcSpectralDifference(denseHRTFdataset_obta,referenceHRTFdataset,[],[],[],true);
titleString = ['OBTA - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 6]);
%saveas(gcf,['SpecDiff_OBTA_N_',num2str(sparseOrder),'.jpg']);

supdeq_calcSpectralDifference(denseHRTFdataset_fdta,referenceHRTFdataset,[],[],[],true);
titleString = ['FDTA - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 6]);
%saveas(gcf,['SpecDiff_FDTA_N_',num2str(sparseOrder),'.jpg']);

supdeq_calcSpectralDifference(denseHRTFdataset_supdeq,referenceHRTFdataset,[],[],[],true);
titleString = ['SUpDEq - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 6]);
%saveas(gcf,['SpecDiff_SUpDEq_N_',num2str(sparseOrder),'.jpg']);

supdeq_calcSpectralDifference(denseHRTFdataset_pc,referenceHRTFdataset,[],[],[],true);
titleString = ['PC - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 6]);
%saveas(gcf,['SpecDiff_PC_N_',num2str(sparseOrder),'.jpg']);

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_pc denseHRTFdataset_supdeq denseHRTFdataset_fdta denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%% Spectral difference in region around the contralateral (only horizontal plane)

testSamplingGrid = [(250:290).',ones(41,1)*90];

supdeq_calcSpectralDifference(denseHRTFdataset_obta,referenceHRTFdataset,testSamplingGrid,[],[],true);
titleString = ['OBTA - Lateral - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 15]);
%saveas(gcf,['SpecDiff_Lateral_OBTA_N_',num2str(sparseOrder),'.jpg']);

supdeq_calcSpectralDifference(denseHRTFdataset_fdta,referenceHRTFdataset,testSamplingGrid,[],[],true);
titleString = ['FDTA - Lateral - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 15]);
%saveas(gcf,['SpecDiff_Lateral_FDTA_N_',num2str(sparseOrder),'.jpg']);

supdeq_calcSpectralDifference(denseHRTFdataset_supdeq,referenceHRTFdataset,testSamplingGrid,[],[],true);
titleString = ['SUpDEq - Lateral - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 15]);
%saveas(gcf,['SpecDiff_Lateral_SUpDEq_N_',num2str(sparseOrder),'.jpg']);

supdeq_calcSpectralDifference(denseHRTFdataset_pc,referenceHRTFdataset,testSamplingGrid,[],[],true);
titleString = ['PC - Lateral - N=',num2str(sparseOrder)];
title(titleString);
ylim([0 15]);
%saveas(gcf,['SpecDiff_Lateral_PC_N_',num2str(sparseOrder),'.jpg']);

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_pc denseHRTFdataset_supdeq denseHRTFdataset_fdta denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%% Get some HRIRs (lateral)

%Define direction
az = 90;
el = 90;

%Get HRIRs
hrir_ref = supdeq_getArbHRIR(referenceHRTFdataset,[az,el]);
hrir_obta = supdeq_getArbHRIR(denseHRTFdataset_obta,[az,el]);
hrir_fdta = supdeq_getArbHRIR(denseHRTFdataset_fdta,[az,el]);
hrir_supdeq = supdeq_getArbHRIR(denseHRTFdataset_supdeq,[az,el]);
hrir_pc = supdeq_getArbHRIR(denseHRTFdataset_pc,[az,el]);

%Plot HRIRs in comparison to reference
supdeq_plotIR(hrir_ref,hrir_obta,[],[],32)
supdeq_plotIR(hrir_ref,hrir_fdta,[],[],32)
supdeq_plotIR(hrir_ref,hrir_supdeq,[],[],32)
supdeq_plotIR(hrir_ref,hrir_pc,[],[],32)

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_pc denseHRTFdataset_supdeq denseHRTFdataset_fdta denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%% Get some HRIRs (contralateral)

%Define direction
az = 270;
el = 90;

%Get HRIRs
hrir_ref = supdeq_getArbHRIR(referenceHRTFdataset,[az,el]);
hrir_obta = supdeq_getArbHRIR(denseHRTFdataset_obta,[az,el]);
hrir_fdta = supdeq_getArbHRIR(denseHRTFdataset_fdta,[az,el]);
hrir_supdeq = supdeq_getArbHRIR(denseHRTFdataset_supdeq,[az,el]);
hrir_pc = supdeq_getArbHRIR(denseHRTFdataset_pc,[az,el]);

%Plot HRIRs in comparison to reference
supdeq_plotIR(hrir_ref,hrir_obta,[],[],32)
supdeq_plotIR(hrir_ref,hrir_fdta,[],[],32)
supdeq_plotIR(hrir_ref,hrir_supdeq,[],[],32)
supdeq_plotIR(hrir_ref,hrir_pc,[],[],32)

%Clear everything besides required variables and final dense HRTF set in SH domain
clearvars -EXCEPT denseHRTFdataset_pc denseHRTFdataset_supdeq denseHRTFdataset_fdta denseHRTFdataset_obta c denseSamplingGrid fs Ndense referenceHRTFdataset rKU100 sparseHRIRdataset sparseHRTFdataset sparseOrder

%% Get circular grids for ITDs/ILDs

%Get circ dataset with 1° in azimuth
circGrid(1:360,1) = 0:359;
circGrid(1:360,2) = 90;

%Get circ grid HRTFs
[circGrid_ref_L,circGrid_ref_R] = supdeq_getArbHRIR(referenceHRTFdataset,circGrid,'DEG',2,'ak'); 
[circGrid_obta_L,circGrid_obta_R] = supdeq_getArbHRIR(denseHRTFdataset_obta,circGrid,'DEG',2,'ak'); 
[circGrid_fdta_L,circGrid_fdta_R] = supdeq_getArbHRIR(denseHRTFdataset_fdta,circGrid,'DEG',2,'ak'); 
[circGrid_supdeq_L,circGrid_supdeq_R] = supdeq_getArbHRIR(denseHRTFdataset_supdeq,circGrid,'DEG',2,'ak'); 
[circGrid_pc_L,circGrid_pc_R] = supdeq_getArbHRIR(denseHRTFdataset_pc,circGrid,'DEG',2,'ak'); 

%% Get ITDs (takes some time)

th = -10; 
for id = 1:size(circGrid,1)
    
    ITD_ref(id) = supdeq_calcITD([circGrid_ref_L(:,id),circGrid_ref_R(:,id)],fs,th);
    ITD_obta(id) = supdeq_calcITD([circGrid_obta_L(:,id),circGrid_obta_R(:,id)],fs,th);
    ITD_fdta(id) = supdeq_calcITD([circGrid_fdta_L(:,id),circGrid_fdta_R(:,id)],fs,th);
    ITD_supdeq(id) = supdeq_calcITD([circGrid_supdeq_L(:,id),circGrid_supdeq_R(:,id)],fs,th);
    ITD_pc(id) = supdeq_calcITD([circGrid_pc_L(:,id),circGrid_pc_R(:,id)],fs,th);
    
end

figure;
lineWidth = 1.2;
plot(ITD_ref,'LineWidth',lineWidth)
hold on;
plot(ITD_obta,'LineWidth',lineWidth); plot(ITD_fdta,'LineWidth',lineWidth); plot(ITD_supdeq,'LineWidth',lineWidth); plot(ITD_pc,'LineWidth',lineWidth);
legend('Ref','OBTA','FDTA','SUpDEq','PA','Location','SouthEast');
xlabel('Azimuth in degrees'); xlim([0 360]);
ylabel('ITD in ms'); 

%% Get broadband ILDs

for id = 1:size(circGrid,1)
    
    ILD_ref(id) = supdeq_calcILD([circGrid_ref_L(:,id),circGrid_ref_R(:,id)]);
    ILD_obta(id) = supdeq_calcILD([circGrid_obta_L(:,id),circGrid_obta_R(:,id)]);
    ILD_fdta(id) = supdeq_calcILD([circGrid_fdta_L(:,id),circGrid_fdta_R(:,id)]);
    ILD_supdeq(id) = supdeq_calcILD([circGrid_supdeq_L(:,id),circGrid_supdeq_R(:,id)]);
    ILD_pc(id) = supdeq_calcILD([circGrid_pc_L(:,id),circGrid_pc_R(:,id)]);
    
end

figure;
lineWidth = 1.2;
plot(ILD_ref,'LineWidth',lineWidth)
hold on;
plot(ILD_obta,'LineWidth',lineWidth); plot(ILD_fdta,'LineWidth',lineWidth); plot(ILD_supdeq,'LineWidth',lineWidth); plot(ILD_pc,'LineWidth',lineWidth);
legend('Ref','OBTA','FDTA','SUpDEq','PC','Location','NorthEast');
xlabel('Azimuth in degrees'); xlim([0 360]);
ylabel('ILD in dB'); 
