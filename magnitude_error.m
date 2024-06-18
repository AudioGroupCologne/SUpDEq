%HRTF_Dataset = load('interpHRTF.mat');

HRTF_Dataset = load('HUTUBS_interp.mat'); % N=3 Sparse Grid


interpHRTF_con = HRTF_Dataset.interpHRTF_con;
interpHRTF_ild = HRTF_Dataset.interpHRTF_ild;
interpHRTF_mca = HRTF_Dataset.interpHRTF_mca;
HRTF_ref = HRTF_Dataset.ref_HRTF;
%ref_HRIR = HRTF_Dataset.ref_HRIR.Data;

% calculate error of the left-ear side
fLowErb = 50;

% con with time alignment
hrirLength = size(interpHRTF_con.HRIR_L,1)/interpHRTF_con.FFToversize;
fs = interpHRTF_con.fs;
[c_con_l,ferb] = AKerbErrorPersistent(interpHRTF_con.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),[fLowErb fs/2],fs);

%Transform to linear values
c_con_l = 10.^(c_con_l/20);

% mca
hrirLength = size(interpHRTF_mca.HRIR_L,1)/interpHRTF_mca.FFToversize;
fs = interpHRTF_mca.fs;
[c_mca_l,~] = AKerbErrorPersistent(interpHRTF_mca.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),...
    [fLowErb fs/2],fs);

%Transform to linear values
c_mca_l = 10.^(c_mca_l/20);

% ild
hrirLength = size(interpHRTF_ild.HRIR_L,1)/interpHRTF_ild.FFToversize;
fs = interpHRTF_ild.fs;

[c_ild_l,~] = AKerbErrorPersistent(interpHRTF_ild.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),...
    [fLowErb fs/2],fs);

%Transform to linear values
c_ild_l = 10.^(c_ild_l/20);

HRTF_Dataset.ref_HRTF.HRIR_L = HRTF_Dataset.ref_HRTF.HRIR_L.';
HRTF_Dataset.ref_HRTF.HRIR_R = HRTF_Dataset.ref_HRTF.HRIR_R.';

%HRTF_Dataset.ref_HRTF.FFToversize = 1;
hrirLength = size(HRTF_Dataset.ref_HRTF.HRIR_L,1)/HRTF_Dataset.ref_HRTF.FFToversize;
fs = HRTF_Dataset.ref_HRTF.f(end)*2;
[c_ref_l,~] = AKerbErrorPersistent(HRTF_Dataset.ref_HRTF.HRIR_L(1:hrirLength,:),AKdirac(hrirLength),...
    [fLowErb fs/2],fs);

%Transform to linear values
c_ref_l = 10.^(c_ref_l/20);

% subject 91
%Get angle for the facing side of the source for ref
azimuthVector_ip_ref = HRTF_Dataset.ref_HRTF.samplingGrid(:,1,1);
iEarsideOrientation_ip_ref = zeros(length(HRTF_Dataset.ref_HRTF.samplingGrid),1);
iEarsideOrientation_ip_ref(azimuthVector_ip_ref == 0 | azimuthVector_ip_ref == 180) = 0; %Center
iEarsideOrientation_ip_ref(azimuthVector_ip_ref > 0 & azimuthVector_ip_ref < 180) = 1; %Left ear
iEarsideOrientation_ip_ref(azimuthVector_ip_ref > 180 & azimuthVector_ip_ref < 360) = -1; %Right ear
averted_ear_ref_idx = find(iEarsideOrientation_ip_ref == 1); %find averted side

%Get angle for the facing side of the source for ref
azimuthVector_ip_mca = HRTF_Dataset.interpHRTF_mca.samplingGrid(:,1,1);
iEarsideOrientation_ip_mca = zeros(length(HRTF_Dataset.interpHRTF_mca.samplingGrid),1);
iEarsideOrientation_ip_mca(azimuthVector_ip_mca == 0 | azimuthVector_ip_mca == 180) = 0; %Center
iEarsideOrientation_ip_mca(azimuthVector_ip_mca > 0 & azimuthVector_ip_mca < 180) = 1; %Left ear
iEarsideOrientation_ip_mca(azimuthVector_ip_mca > 180 & azimuthVector_ip_mca < 360) = -1; %Right ear
averted_ear_mca_idx = find(iEarsideOrientation_ip_mca == 1); %find averted side


% overall error mca 
%dG_mca = abs(10*log10(mean(c_mca_l,"all")) / mean(c_ref_l,"all"));
dG_mca = abs(10*log10(mean(c_mca_l(:,averted_ear_mca_idx),"all")) / mean(c_ref_l(:,averted_ear_ref_idx),"all"));
disp(['ΔG mca: ',num2str(dG_mca), ' dB mit N = 3 '])

%dG_ild = abs(10*log10(mean(c_ild_l,"all")) / mean(c_ref_l,"all"));
dG_ild = abs(10*log10(mean(c_ild_l(:,averted_ear_mca_idx),"all")) / mean(c_ref_l(:,averted_ear_ref_idx),"all"));
disp(['ΔG ild: ',num2str(dG_ild), ' dB mit N = 3 '])


%average over frequency and (subjects) deltaG(omega)
%dG_avg_freq_subj = abs(10*log10(mean(c_mca_l,1)) ./ mean(c_ref_l,1));

%average over direction deltaG(f_c,s)
%dG_avg_dir = abs(10*log10(mean(c_mca_l,2)) ./ mean(c_ref_l,2));

%disp(['ΔG: ',num2str(dG), ' dB mit N = ',num2str(ref_HRTF.N)])
%disp(['ΔG: ',num2str(dG), ' dB mit N = 3 ')

% save the magnitude error data
%save('magnitude_error','dG','dG_avg_freq_subj','dG_avg_dir')

