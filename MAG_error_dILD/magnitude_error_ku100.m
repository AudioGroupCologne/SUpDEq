%% Calculate total Magnitude error ΔG and over directions ΔG(Ω), frequencies ΔG(f)

clear all;

N = input('Enter the SH-Order the interpolated Sparse Grid (1-5): ');
HRTF_Dataset = load([fullfile('mat_data','ku100_interp_N'),num2str(N),'.mat']);

%Get interpolated and Reference dataset
%interpHRTF_con = HRTF_Dataset.interpHRTF_con;
interpHRTF_ild = HRTF_Dataset.interpHRTF_ild;
interpHRTF_mca = HRTF_Dataset.interpHRTF_mca;
HRTF_ref = HRTF_Dataset.ref_HRTF;

%Smooth interpolated HRTF's (left channel)in Auditory Bands 
%to calaculate the magnitude error

%CTA - dataset
%hrirLength = size(interpHRTF_con.HRIR_L,1)/interpHRTF_con.FFToversize;
%fs = interpHRTF_con.fs;l
%[c_con_l,ferb] = AKerbErrorPersistent(interpHRTF_con.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),[fLowErb fs/2],fs);
%Transform to linear values
%c_con_l = 10.^(c_con_l/20);

%MCA - dataset
fLowErb = 50;
hrirLength = size(interpHRTF_mca.HRIR_L,1)/interpHRTF_mca.FFToversize;
fs = interpHRTF_mca.fs;
[c_mca_l,~] = AKerbErrorPersistent(interpHRTF_mca.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),[fLowErb fs/2],fs);
%Transform to linear values
c_mca_l = 10.^(c_mca_l/20);

%ILD - dataset
hrirLength = size(interpHRTF_ild.HRIR_L,1)/interpHRTF_ild.FFToversize;
fs = interpHRTF_ild.fs;
[c_ild_l,~] = AKerbErrorPersistent(interpHRTF_ild.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),[fLowErb fs/2],fs);
%Transform to linear values
c_ild_l = 10.^(c_ild_l/20);

%Reference - dataset
hrirLength = size(HRTF_ref.HRIR_L,1)/HRTF_ref.FFToversize;
fs = HRTF_ref.f(end)*2;
[c_ref_l,~] = AKerbErrorPersistent(HRTF_ref.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),[fLowErb fs/2],fs);
%Transform to linear values
c_ref_l = 10.^(c_ref_l/20);

%Calculate total magnitude error ΔG of the left channel 
dG_mca = mean(abs(10*log10(c_mca_l ./ c_ref_l)),'all');
dG_ild = mean(abs(10*log10(c_ild_l ./ c_ref_l)),'all');

%Average over frequencies --> ΔG(Ω)
dG_avg_mca_freq = mean(abs(10*log10(c_mca_l ./ c_ref_l)),1);
dG_avg_ild_freq = mean(abs(10*log10(c_ild_l ./ c_ref_l)),1);

%Average over direction --> ΔG(f)
dG_avg_mca_dir = mean(abs(10*log10(c_mca_l ./ c_ref_l)),2);
dG_avg_ild_dir = mean(abs(10*log10(c_ild_l ./ c_ref_l)),2);

%Display the total magnitude error of the HRTF's (left channel) with only
%applied mca-correction and once with mca-correction and ild-correction
disp(['ΔG mca: ',num2str(dG_mca), ' dB with N = ',num2str(N)])
disp(['ΔG ild: ',num2str(dG_ild), ' dB with N = ',num2str(N)])


%% Save Magnitude_errors in .mat files for further processing
save([fullfile('mat_data','ku100_magnitude_error_N'),num2str(N)],'dG_mca','dG_avg_mca_freq','dG_avg_mca_dir','dG_ild','dG_avg_ild_freq','dG_avg_ild_dir')
