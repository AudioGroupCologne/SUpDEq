function [] = plot_llado2022(y_est_dir,y_est_uncertainty,angle_id,y_test)
%PLOT_LLADO2022 plots the figures from llado et al. (2022)
%   Usage: plot_llado2022(y_est_dir,y_est_uncertainty,angle_id,y_test);
%
%   Input parameters:
%     y_est_dir          : Estimated data for perceived direction
%     y_est_uncertainty  : Estimated data for perceived uncertainty
%     angle_id           : Vector of lateral angles of the test subset
%
%   PLOT_LLADO2022 plots the Figures from llado2022.
%
%   Optional input parameters:
%
%     'y_test'           labels from subjective test
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_llado2022.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB M - Communication Systems
%   #Author: Pedro Llado (2022)
%   #Author: Petteri Hyvärinen (2022)
%   #Author: Ville Pulkki (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%%  DEFAULT OPTIONAL INPUTS
%%
figure;
subplot(1,2,1)
if nargin > 3
    plot(angle_id,y_test(:,1),'k','LineWidth',2);
    hold on;
end
plot(angle_id,y_est_dir,'k','LineStyle','--','LineWidth',2);
title('Perceived direction');
ylim([min(angle_id)-10 max(angle_id)+10])
xlim([min(angle_id)-10 max(angle_id)+10])
ax = gca();
grid(ax, 'on') 
set(ax,'XTick',angle_id)
ylabel('Perceived direction (°)')
xlabel('Sound source direction (°)')

ax.XTick = [angle_id];
ax.XTickLabel = [angle_id];


ax.YTick = [angle_id];
ax.YTickLabel = [angle_id(end:-1:1)];
set(gca, 'YDir','reverse')
if nargin > 3
    legend('Subjective','NN-Estimated','Location','northwest');
else
    legend('NN-Estimated','Location','northwest');
end



%%
%figure;
subplot(1,2,2)
if nargin >3
    plot(angle_id,y_test(:,2),'k','LineWidth',2);
    hold on;
end
plot(angle_id,y_est_uncertainty,'k','LineStyle','--','LineWidth',2);
title('Position uncertainty');
ylim([0 10])
xlim([min(angle_id)-10 max(angle_id)+10])
ax = gca();
grid(ax, 'on') 
set(ax,'XTick',angle_id)
ylabel('Position uncertainty (°)')
xlabel('Sound source direction (°)')
ax.XTick = [angle_id];
ax.XTickLabel = angle_id;

if nargin > 3
    legend('Subjective','NN-Estimated','Location','northwest');
else
    legend('NN-Estimated','Location','southwest');
end
end


