%% (8) Compare HRTF of the Reference set with the interpolated HRTF of ^H and the magnitude-corrected interpolated HRTF of ^H_C with the applied ILD-Filter

% ref = SHcoeffficiets_simulated;
% samplingGrid = sgD;
% mode = 'DEG';
% channel = 2;
% transformCore = 'ak';
% [denseHRTFdataset.HRTF_L, denseHRTFdataset.HRTF_R] = supdeq_getArbHRTF_HUTUBS(ref,samplingGrid,mode,channel,transformCore);
% 
% NmaxSparse= Ns;
% FFToversize = 1; % ----> figure out the right fftoversize ...
% denseGrid = sgD;
% 
% %Fill struct with additional info
% denseHRTFdataset.f = ref.HRIR.f;
% denseHRTFdataset.fs = ref.HRIR.f(end)*2;
% denseHRTFdataset.Nmax = NmaxSparse;
% denseHRTFdataset.FFToversize = FFToversize;
% denseHRTFdataset.samplingGrid = sparseGrid;
% if isfield(ref,'ReceiverPosition')
%     denseHRTFdataset.sourceDistance = ref.ReceiverPosition;
% end
% 
% ref_HRTF = denseHRTFdataset;

%Get reference data
FFToversize_ref = 1; %Set FFT-Blocksize
sparseHRIRdataset_SOFA = HRIRdataset_simulated_SOFA;
%Convert reference  hrir to hrtf
ref_HRTF = supdeq_sofa2hrtf(sparseHRIRdataset_SOFA, Nd,[],FFToversize);

% add SH Order to struct for further processing
ref_HRTF.N = Nd;

%Get SampelingGrid of the ref and interpolated data
ref_HRTF_sampgrid = ref_HRTF.samplingGrid(:,1:2);
intpHRTF_con_sampgrid = interpHRTF_con.samplingGrid(:,1:2);
intpHRTF_mca_sampgrid = interpHRTF_mca.samplingGrid(:,1:2);
intpHRTF_ild_sampgrid = interpHRTF_ild.samplingGrid(:,1:2);

% set coordinate_system --> grids are generate with [180° 90° 0°]
%coordinate_system = input('Enter 0 for elevation: [90° 0° -90°] and enter 1 for elevation: [180° 90° 0°]: '); % N:180° S:0°

%default angle:
% azimuth: [0, 90, 180, 270]
% elevation: [90, 0 -90]

% --generate grid with supdeq_lebedev() or supdeq_fliege() -> always [180°90°0°]--

% --> enter 1 to adapt the angles to [90° 0° -90°]
%coordinate_system = input('Enter 0 for elevation: [90° 0° -90°] 
% and enter 1 for elevation: [180° 90° 0°]: '); % N:180° S:0°


%if coordinate_system
    %Adjust Elevation angle according to Magnitude-Corrected and Time-Aligned
    %Interpolation of Head-Related Transfer Functions Paper if
    %Dataset-coordinates are differ
    ref_HRTF_sampgrid(:,2) = ref_HRTF.samplingGrid(:,2)-90;
    intpHRTF_con_sampgrid(:,2) = interpHRTF_con.samplingGrid(:,2)-90;
    intpHRTF_mca_sampgrid(:,2) = interpHRTF_mca.samplingGrid(:,2)-90;
    intpHRTF_ild_sampgrid(:,2) = interpHRTF_ild.samplingGrid(:,2)-90;

%end

%Set azimuth and elevation angle
azimuth = input('Enter azimuth angle: ');
elevation = input('Enter elevation angle: ');

%Set Target angle [azimuth elevation]
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

%Plot the reference, conv. time aligned and mca + ILD-Filter HRTF's at
% desired coordinates to compare

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

xlabel('Frequency in Hz','LineWidth',14);
ylabel('Magnitude in dB','LineWidth',14);
title('HUTUBS Subject 91')
grid on;


%xlim([interpHRTF_mca.f(1) interpHRTF_mca.f(end)])
%xlim([interpHRTF_ild.f(1) interpHRTF_ild.f(end)])

legend([p1,p3,p5,p7],'$H_{\mathrm{R}}$','$\widehat{H}$','$\widehat{H}_{\mathrm{ILD}}$','$\widehat{H}_{\mathrm{C}}$', ...
    'Interpreter','latex','FontSize',20,'Location','southwest');

xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');

str = sprintf(' = (%.2f, %.2f)',[intpHRTF_mca_sampgrid(idx_intpHRTF_mca,1),intpHRTF_mca_sampgrid(idx_intpHRTF_mca,2)]);
combinedStr = strcat('$\Omega$ ',str);
text(xlim(1)+250, ylim(2)-3,combinedStr ,'Interpreter','latex','FontSize',20,'FontWeight','bold')

%combinedStr = ['N = ',num2str(Nd)];
%text(xlim(1)+250, ylim(2)-3,combinedStr ,'Interpreter','latex','FontSize',20,'FontWeight','bold')

%Set R and L text
text(interpHRTF_mca.f(idx_max_mc), max_interpHRTF_mca + 1.5,upper_lab,'FontSize',24,'FontWeight','bold')
text(interpHRTF_mca.f(idx_max_mc), value_other_channel - 1.5,lower_lab,'FontSize',24,'FontWeight','bold')
