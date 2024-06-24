%% Plot magnitude errors ΔG(Ω) over directions Ω

% Request user input for SH-Order
N = input('Enter the SH-Order the interpolated Sparse Grid (1-5): ');

%Load data
HRTF_Dataset = load([fullfile('mat_data','ku100_interp_N'), num2str(N), '.mat']);
mag_error = load([fullfile('mat_data','ku100_magnitude_error_N'), num2str(N), '.mat']);

%Extract data
X_data = HRTF_Dataset.interpHRTF_ild.samplingGrid(:, 1); % Azimuth
Y_data = HRTF_Dataset.interpHRTF_ild.samplingGrid(:, 2); % Elevation
data_mca = mag_error.dG_avg_mca_freq; 
data_ild = mag_error.dG_avg_ild_freq; % Magnitude error data

%Create grid
[x, ~, xi] = unique(X_data);
[y, ~, yi] = unique(Y_data);
subs = [xi yi]; %Use the indices instead of values as coordinates


%% Plot magnitude errors MCA HRTF

%Fill the grid with data and handle missing values
M = accumarray(subs, data_mca, [], [], NaN); %Fill missing with NaN
[X, Y] = ndgrid(x, y); %Generate the grid surf() expects

%Interpolate missing values
Z = fillmissing(M, 'linear', 'EndValues', 'nearest'); % Alternatives: 'nearest' or 'spline'

figure(1)
surf(X, Y, Z)
view(0, 90)
xlabel('Azimuth in degree', 'FontSize', 14) % Increased font size
ylabel('Elevation in degree', 'FontSize', 14) % Increased font size
xticks(0:30:360);
yticks(-90:30:90);
xlim([0 361])

%Set title
title(['ku100 with only Mag. corr. N = ', num2str(N),',   \DeltaG=', num2str(round(mag_error.dG_mca,2)), ' dB'],FontSize=14)

%Colorbar and colormap adjustments
colorbar
colormap('hot') %Adjust colormap to match the color scheme
clim([0, 5]);
shading interp %Interpolated shading for smooth color transitions
set(gca, 'YDir', 'normal') %Ensure y-axis is in normal direction
axis tight %Fit the axes to the data

%Additional adjustments for better visualization
set(gca, 'FontSize', 12) % Increase font size for tick labels



%% Plot magnitude errors ILD HRTF

%Fill the grid with data and handle missing values
M = accumarray(subs, data_ild, [], [], NaN); %Fill missing with NaN
[X, Y] = ndgrid(x, y); %Generate the grid surf() expects

%Interpolate missing values
Z = fillmissing(M, 'linear', 'EndValues', 'nearest'); % Alternatives: 'nearest' or 'spline'


%Plot ILD HRTF
figure(2)
surf(X, Y, Z)
view(0, 90)
xlabel('Azimuth in degree', 'FontSize', 14) % Increased font size
ylabel('Elevation in degree', 'FontSize', 14) % Increased font size
xticks(0:30:360);
yticks(-90:30:90);
xlim([0 361])

title(['ku100 with Mag. corr. and ILD Comp N = ', num2str(N),',   \DeltaG=', num2str(round(mag_error.dG_ild,2)), ' dB'],FontSize=14)

%Colorbar and colormap adjustments
colorbar
colormap('hot') % Adjust colormap to match the color scheme
clim([0, 5]);
shading interp %Interpolated shading for smooth color transitions
set(gca, 'YDir', 'normal') %Ensure y-axis is in normal direction
axis tight %Fit the axes to the data

%Additional adjustments for better visualization
set(gca, 'FontSize', 12) % Increase font size for tick labels

