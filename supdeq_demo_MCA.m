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

clear variables;
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
% 
%Plot frontal left-ear HRIR (Az = 0, Col = 90) of conventional, MCA and ILD-Filter
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

% Left ear
%[erb_sh,fc_erb] = AKerbError(interpHRTF_sh.HRIR_L,refHRTF.HRIR_L, [50 fs/2], fs);
%erb_con = AKerbError(interpHRTF_con.HRIR_L,refHRTF.HRIR_L, [50 fs/2], fs);
%erb_mca = AKerbError(interpHRTF_mca.HRIR_L,refHRTF.HRIR_L, [50 fs/2], fs);

% Right ear
[erb_sh,fc_erb] = AKerbError(interpHRTF_sh.HRIR_R,refHRTF.HRIR_R, [50 fs/2], fs);
erb_con = AKerbError(interpHRTF_con.HRIR_R,refHRTF.HRIR_R, [50 fs/2], fs);
erb_mca = AKerbError(interpHRTF_mca.HRIR_R,refHRTF.HRIR_R, [50 fs/2], fs);

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


%% (8) Compare HRTF of the Reference set with the interpolated HRTF of ^H and the magnitude-corrected interpolated HRTF of ^H_C with the applied ILD-Filter

target = [60, 60]; % target values azimuth; elevation angle

refHRTF_sampgrid = refHRTF.samplingGrid(:,1:2);
intpHRTF_con_sampgrid = interpHRTF_con.samplingGrid(:,1:2);
intpHRTF_mca_sampgrid = interpHRTF_mca.samplingGrid(:,1:2);

% find the nearest coordinates next to the target values
diff_refHRTF = abs(refHRTF_sampgrid - target);
diff_intpHRTF_con = abs(intpHRTF_con_sampgrid - target);
diff_intpHRTF_mca = abs(intpHRTF_mca_sampgrid - target);

sum_diff_refHRTF = sum(diff_refHRTF,2);
sum_diff_intpHRTF_con = sum(diff_intpHRTF_con,2);
sum_diff_intpHRTF_mca = sum(diff_intpHRTF_mca,2);

[~, idx_refHRTF] = min(sum_diff_refHRTF);
[~, idx_intpHRTF_con] = min(sum_diff_intpHRTF_con);
[~, idx_intpHRTF_mca] = min(sum_diff_intpHRTF_mca);

% recognize the facing side (channel) to the source
con_HRTF_L = 20*log10(abs(interpHRTF_con.HRTF_L(idx_intpHRTF_con,:)));
con_HRTF_R = 20*log10(abs(interpHRTF_con.HRTF_R(idx_intpHRTF_con,:)));

% determine the energetic sum
L_sum_HRTF_L = 10 * log10(sum(10.^(con_HRTF_L / 10))); 
L_sum_HRTF_R = 10 * log10(sum(10.^(con_HRTF_R / 10)));

if L_sum_HRTF_L > L_sum_HRTF_R

    %[max_ref,idx_refHRT] = max(20*log10(abs(refHRTF.HRTF_L(idx_refHRTF,:))));
    %[max_con,idx_max_con] = max(20*log10(abs(interpHRTF_con.HRTF_L(idx_refHRTF,:))));
    [max_mc,idx_max_mc] = max(20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,:))));
    
    max_interpHRTF_mca = 20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,idx_max_mc)));
    value_other_channel = 20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,idx_max_mc)));
    upper_lab = 'L';
    lower_lab = 'R';

elseif L_sum_HRTF_L < L_sum_HRTF_R

   %[max_ref,idx_refHRT] = max(20*log10(abs(refHRTF.HRTF_R(idx_refHRTF,:))));
   %[max_con,idx_max_con] = max(20*log10(abs(interpHRTF_con.HRTF_R(idx_refHRTF,:))));
   [max_mc,idx_max_mc] = max(20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:))));
   
   max_interpHRTF_mca = 20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,idx_max_mc)));
   value_other_channel = 20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,idx_max_mc)));
   upper_lab = 'R';
   lower_lab = 'L';

end
% Plot the reference, conv. time aligned and mca + ILD-Filter HRTF's at
% desired coordinates to compare

AKf(18,9);
% Reference HRTF
p1 = semilogx(refHRTF.f,20*log10(abs(refHRTF.HRTF_L(idx_refHRTF,:))) ...
    ,'LineWidth',2.0,'Color',[0,0,66]/255); %,'alpha', 0.8);
hold on;
p2 = semilogx(refHRTF.f,20*log10(abs(refHRTF.HRTF_R(idx_refHRTF,:))) ...
    ,'LineWidth',2.0,'Color',[0,0,66]/255);  %,'alpha', 0.8);

% conventional HRTF 
p3 = semilogx(interpHRTF_con.f,20*log10(abs(interpHRTF_con.HRTF_L(idx_intpHRTF_con,:))), ... 
    'LineWidth',2.0,'Color',[160,160,160]/255);
p4 = semilogx(interpHRTF_con.f,20*log10(abs(interpHRTF_con.HRTF_R(idx_intpHRTF_con,:))), ... 
    'LineWidth',2.0,'Color',[160,160,160]/255);

% mca and ILD corrected HRTF
p5 = semilogx(interpHRTF_mca.f,20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,:))) ... 
    ,'LineWidth',2.0,'Color',[255,102,102]/255);
p6 = semilogx(interpHRTF_mca.f,20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:))) ... 
    ,'LineWidth',2.0,'Color',[255,102,102]/255);

xlabel('Frequency in Hz','LineWidth',14);
ylabel('Magnitude in dB','LineWidth',14);
grid on;
xlim([100 20000])

legend([p1,p3,p5],'$H_{\mathrm{R}}$','$\widehat{H}$','$\widehat{H}_{\mathrm{C}}$', ...
    'Interpreter','latex','FontSize',20,'Location','southwest');

xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');

str = sprintf(' = (%.2f, %.2f)',[intpHRTF_mca_sampgrid(idx_intpHRTF_mca,1),intpHRTF_mca_sampgrid(idx_intpHRTF_mca,2)]);
combinedStr = strcat('$\Omega$ ',str);
text(xlim(1)+10, ylim(2)-3,combinedStr ,'Interpreter','latex','FontSize',20,'FontWeight','bold')

% set R and L text
text(interpHRTF_mca.f(idx_max_mc), max_interpHRTF_mca + 1.5,upper_lab,'FontSize',24,'FontWeight','bold')
text(interpHRTF_mca.f(idx_max_mc), value_other_channel - 1.5,lower_lab,'FontSize',24,'FontWeight','bold')



