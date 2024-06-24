%% (9) - Plot and compare Reference HRTF's with CTA-HRTF's, MCA-HRFT's, ILD-HRTF'S
%Use this script, to generate the datasets by supdey_demo_MCA.m
clear all;

%Request user input for SH-Order
Ns = input('Enter the SH-Order the interpolated Sparse Grid (1-5): ');

%Load data
HRTF_Dataset = load([fullfile('mat_data','ku100_interp_N'), num2str(Ns), '.mat']);

%Get the datasets
interpHRTF_con = HRTF_Dataset.interpHRTF_con;
interpHRTF_mca = HRTF_Dataset.interpHRTF_mca;
interpHRTF_ild = HRTF_Dataset.interpHRTF_ild;
ref_HRTF = HRTF_Dataset.ref_HRTF;

%Get SampelingGrid of the Refernce and interpolated data
ref_HRTF_sampgrid = round(ref_HRTF.samplingGrid(:,1:2));
intpHRTF_con_sampgrid = round(interpHRTF_con.samplingGrid(:,1:2));
intpHRTF_mca_sampgrid = round(interpHRTF_mca.samplingGrid(:,1:2));
intpHRTF_ild_sampgrid = round(interpHRTF_ild.samplingGrid(:,1:2));

%Set azimuth and elevation angle
azimuth = input('Enter azimuth angle: ');
elevation = input('Enter elevation angle: ');

%Set target angle [azimuth elevation]
target = [azimuth, elevation];

%Get angle for the facing side of the source
azimuthVector_ip = ref_HRTF.samplingGrid(:,1,1);
iEarsideOrientation_ip = zeros(length(ref_HRTF.samplingGrid),1);
iEarsideOrientation_ip(azimuthVector_ip == 0 | azimuthVector_ip == 180) = 0; %Center
iEarsideOrientation_ip(azimuthVector_ip > 0 & azimuthVector_ip < 180) = 1; %Left ear
iEarsideOrientation_ip(azimuthVector_ip > 180 & azimuthVector_ip < 360) = -1; %Right ear

%Find the nearest coordinates next to the target values
%(Should always be the same idx only for test purposes)
diff_refHRTF = sum(abs(ref_HRTF_sampgrid - target),2);
diff_intpHRTF_con = sum(abs(intpHRTF_con_sampgrid - target),2);
diff_intpHRTF_mca = sum(abs(intpHRTF_mca_sampgrid - target),2);
diff_intpHRTF_ild = sum(abs(intpHRTF_ild_sampgrid - target),2);

%Get index of the target angle
[~, idx_refHRTF] = min(diff_refHRTF);
[~, idx_intpHRTF_con] = min(diff_intpHRTF_con);
[~, idx_intpHRTF_mca] = min(diff_intpHRTF_mca);
[~, idx_intpHRTF_ild] = min(diff_intpHRTF_ild);

%Set L & R marker
if iEarsideOrientation_ip(idx_refHRTF) == 1 %Left ear

    [max_mc,idx_max_mc] = max(20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,:))));
    max_interpHRTF_mca = 20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,idx_max_mc)));
    value_other_channel = 20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,idx_max_mc)));
    upper_lab = 'L';
    lower_lab = 'R';

elseif iEarsideOrientation_ip(idx_refHRTF) == -1 %Right ear

    [max_mc,idx_max_mc] = max(20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:))));
    max_interpHRTF_mca = 20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,idx_max_mc)));
    value_other_channel = 20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,idx_max_mc)));
    upper_lab = 'R';
    lower_lab = 'L';

elseif iEarsideOrientation_ip(idx_refHRTF) == 0 %Center --> test with center == L
    [max_mc,idx_max_mc] = max(20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:))));
    max_interpHRTF_mca = 20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,idx_max_mc)));
    value_other_channel = 20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,idx_max_mc)));
    upper_lab = 'R';
    lower_lab = 'L';
end

%Plot HRTF's
AKf(18,9);

%Reference HRTF
p1 = semilogx(ref_HRTF.f,20*log10(abs(ref_HRTF.HRTF_L(idx_refHRTF,:))) ...
    ,'LineWidth',2.0,'Color',[0,0,66]/255); %,'alpha', 0.8);
hold on;
p2 = semilogx(ref_HRTF.f,20*log10(abs(ref_HRTF.HRTF_R(idx_refHRTF,:))) ...
    ,'LineWidth',2.0,'Color',[0,0,66]/255);  %,'alpha', 0.8);

%Conventional HRTF
p3 = semilogx(interpHRTF_con.f,20*log10(abs(interpHRTF_con.HRTF_L(idx_intpHRTF_con,:))), ...
    'LineWidth',2.0,'Color',[160,160,160]/255);
p4 = semilogx(interpHRTF_con.f,20*log10(abs(interpHRTF_con.HRTF_R(idx_intpHRTF_con,:))), ...
    'LineWidth',2.0,'Color',[160,160,160]/255);

%Mca HRTFT + ILD-filter on the averted side
p5 = semilogx(interpHRTF_ild.f,20*log10(abs(interpHRTF_ild.HRTF_L(idx_intpHRTF_ild,:))) ...
    ,'LineWidth',2.0,'Color',[0,255,0]/255);
p6 = semilogx(interpHRTF_ild.f,20*log10(abs(interpHRTF_ild.HRTF_R(idx_intpHRTF_ild,:))) ...
    ,'LineWidth',2.0,'Color',[0,255,0]/255);

%Mca HRTFT
p7 = semilogx(interpHRTF_mca.f,20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,:))) ...
    ,'LineWidth',2.0,'Color',[255,102,102]/255);
p8 = semilogx(interpHRTF_mca.f,20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:))) ...
    ,'LineWidth',2.0,'Color',[255,102,102]/255);

%Set xlabel,ylabel, title ...
xlabel('Frequency in Hz','LineWidth',14);
ylabel('Magnitude in dB','LineWidth',14);

title('KU100 HRTF','LineWidth',16)
grid on;
clear xlim;
clear ylim;

xlim([interpHRTF_ild.f(1) interpHRTF_ild.f(end)]);

if strcmp(upper_lab,'L')
    ylim([min(20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:)))) ...
        + 0.2 * min(20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:)))) ...
        max(20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,:)))) + 5]);
else
    ylim([min(20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,:)))) ...
        + 0.2 * min(20*log10(abs(interpHRTF_mca.HRTF_L(idx_intpHRTF_mca,:)))) ...
        max(20*log10(abs(interpHRTF_mca.HRTF_R(idx_intpHRTF_mca,:)))) + 5]);
end

xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');

%Set legend
legend([p1,p3,p5,p7],'$H_{\mathrm{R}}$','$\widehat{H}$','$\widehat{H}_{\mathrm{ILD}}$','$\widehat{H}_{\mathrm{C}}$', ...
    'Interpreter','latex','FontSize',20,'Location','southwest');

%Extract and round values
values = [intpHRTF_mca_sampgrid(idx_intpHRTF_mca,1), intpHRTF_mca_sampgrid(idx_intpHRTF_mca,2)];
rounded_values = round(values, 0);

%Combine strings
combinedStr = strcat('$\Omega =$','(',num2str(rounded_values(1)),',',num2str(rounded_values(2)),')');
text(120, ylim(2)-3,combinedStr ,'Interpreter','latex','FontSize',20,'FontWeight','bold')
combinedStr = ['N = ',num2str(Ns)];
text(200, ylim(2)-3,combinedStr ,'Interpreter','latex','FontSize',20,'FontWeight','bold')

%Set R and L
text(interpHRTF_mca.f(idx_max_mc), max_interpHRTF_mca + 1.5,upper_lab,'FontSize',24,'FontWeight','bold')
text(interpHRTF_mca.f(idx_max_mc), value_other_channel - 1.5,lower_lab,'FontSize',24,'FontWeight','bold')