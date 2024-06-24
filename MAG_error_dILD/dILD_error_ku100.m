%% Calculate total ILD error ΔILD and ILD errors over directions Ω ΔILD(Ω)

clear all;

N = input('Enter the SH-Order the interpolated Sparse Grid (1-5): ');
HRTF_Dataset = load([fullfile('mat_data','ku100_interp_N'),num2str(N),'.mat']); 

%Get interpolated dataset
%interpHRTF_con = HRTF_Dataset.interpHRTF_con;
interpHRTF_ild = HRTF_Dataset.interpHRTF_ild;
interpHRTF_mca = HRTF_Dataset.interpHRTF_mca;
HRTF_ref = HRTF_Dataset.ref_HRTF;

%To calculate the absolute ILD errors only directions Ω in the horizontal plane are needed 
%find only grid points with elevation angle = 0° if desired tolerance of
%elevation angle can be increased --> more grid points will be found
%azimuth: 0° <= phi <= 359, elevation: 0°

%Get the idxes of the gridpoints in the horizontal plane
tolerance = 0; 
elevationVector_ip = interpHRTF_mca.samplingGrid(:,2);
elevation_side_idx = find(elevationVector_ip >= (0 - tolerance) & elevationVector_ip <= (0 + tolerance));
ILD_grid = interpHRTF_mca.samplingGrid(elevation_side_idx,:);

%Calculate the absolute ILD_errors in dependence of the angle in horizontal
%plane ΔILD(Ω) of HRIR with magnitude correction and HRTF's with both
%magnitude correction and ild correction
d_ILD_mca_ang = abs(10*log10(sum(interpHRTF_mca.HRIR_L(:,elevation_side_idx).^2,1) ./ sum(interpHRTF_mca.HRIR_R(:,elevation_side_idx).^2,1)) ...
       - 10*log10(sum(HRTF_ref.HRIR_L(:,elevation_side_idx).^2,1) ./ sum(HRTF_ref.HRIR_R(:,elevation_side_idx).^2,1)));

d_ILD_ild_ang = abs(10*log10(sum(interpHRTF_ild.HRIR_L(:,elevation_side_idx).^2,1) ./ sum(interpHRTF_ild.HRIR_R(:,elevation_side_idx).^2,1)) ...
    - 10*log10(sum(HRTF_ref.HRIR_L(:,elevation_side_idx).^2,1) ./ sum(HRTF_ref.HRIR_R(:,elevation_side_idx ).^2,1)));

%Calculate the total absolute ILD_errors averaged over directions Ω
d_ILD_mca_ges = mean(d_ILD_mca_ang);
d_ILD_ild_ges = mean(d_ILD_ild_ang);

%Display the total absolute ILD_errors averaged over directions Ω
disp(['ΔILD mca: ',num2str(d_ILD_mca_ges), ' dB with N = ',num2str(N)])
disp(['ΔILD ild: ',num2str(d_ILD_ild_ges), ' dB with N = ',num2str(N)])


%% Save ILD_errors in .mat file for further processing/calculations
save([fullfile('mat_data','ku100_dILD_error_N'),num2str(N)],'d_ILD_mca_ges','d_ILD_mca_ang','d_ILD_ild_ges','d_ILD_ild_ang','ILD_grid')
