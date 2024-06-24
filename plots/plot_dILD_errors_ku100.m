%% PLOT dILD errors over directions Ω

clear all;

%Request user input for SH-Order
N = input('Enter the SH-Order the interpolated Sparse Grid (1-5): ');

%Load data
%HRTF_Dataset = load([fullfile('mat_data','ku100_interp_N'),num2str(N),'.mat']);
ILD_error = load([fullfile('mat_data','ku100_dILD_error_N'), num2str(N), '.mat']);

%adapt the azimuth grid dor visulation porpuses
azimuth_angle =  ILD_error.ILD_grid(:,1) -180;

%% Plot dILD errors MCA-HRTF's

figure(1)
plot(fftshift(circshift(azimuth_angle,int32(length(ILD_error.ILD_grid(:,1))/2)))...
    ,fftshift(circshift(ILD_error.d_ILD_mca_ang.',int32(length(ILD_error.d_ILD_mca_ang.')/2))),linewidth=2)

hold on

% Plot JND line
plot(fftshift(circshift(azimuth_angle,int32(length(ILD_error.ILD_grid(:,1))/2))), 0.93 * ones(size(azimuth_angle)), 'r--',linewidth=1.2);

grid on

% Get current axes limits
xlim([-180.5 180.5])
ylim([0 10.5])
x_limits = xlim;
y_limits = ylim;

xlabel('Azimuth in degree','FontSize',14)
ylabel('ΔILD (Ω) in dB','FontSize',14)

title(['ku100 only with Mag. corr. ILD errors N = ',num2str(N)],'FontSize',14);

% Set the xtick values
xticks(-180:30:180);

legend('ku100','JND',fontsize=14)
% Calculate position for the text (near the upper left corner)
text_x = x_limits(1) + (x_limits(2) - x_limits(1)) * 0.05; % 5% from the left side
text_y = y_limits(2) * 0.95; % 95% of the y-axis limit

% Add the text label
text(text_x, text_y,['N = ',num2str(N), ',  ΔILD = ', num2str(round(ILD_error.d_ILD_mca_ges,2)),' dB'], 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'left');

%% Plot dILD errors ILD-HRTF's

figure(2)
plot(fftshift(circshift(azimuth_angle,int32(length(ILD_error.ILD_grid(:,1))/2)))...
    ,fftshift(circshift(ILD_error.d_ILD_ild_ang.',int32(length(ILD_error.d_ILD_ild_ang.')/2))),linewidth=2)

hold on

% Plot JND line
plot(fftshift(circshift(azimuth_angle,int32(length(ILD_error.ILD_grid(:,1))/2))), 0.93 * ones(size(azimuth_angle)), 'r--',linewidth=1.2);

grid on

% Get current axes limits
xlim([-180.5 180.5])
ylim([0 10.5])
x_limits = xlim;
y_limits = ylim;
title(['ku100 delta ILD mca N',num2str(N)]);

xlabel('Azimuth in degree','FontSize',14)
ylabel('ΔILD (Ω) in dB','FontSize',14)

title(['ku100 with Mag. corr. and ILD corr. ILD errors N = ',num2str(N)],'FontSize',14);

% Set the xtick values
xticks_vector = -180:30:180;
xticks(xticks_vector);

legend('ku100','JND',fontsize=14)
% Calculate position for the text (near the upper left corner)
text_x = x_limits(1) + (x_limits(2) - x_limits(1)) * 0.05; % 5% from the left side
text_y = y_limits(2) * 0.95; % 95% of the y-axis limit

% Add the text label
text(text_x, text_y,['N = ',num2str(N), ',  ΔILD = ', num2str(round(ILD_error.d_ILD_ild_ges,2)),' dB'], 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'left');
