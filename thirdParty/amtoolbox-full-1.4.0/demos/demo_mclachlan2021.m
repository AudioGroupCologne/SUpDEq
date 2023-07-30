%DEMO_MCLACHLAN2021 Demo for dynamic full sphere localization model based on Reijniers et al. (2014)
%
%   DEMO_MCLACHLAN2021(flag) demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   on the full sphere using the localization model from Reijniers et al. (2014).
%
%   Figure 1: Baseline prediction: averaged polar and lateral accuracy
%
% 
%   This demo computes the baseline prediction (localizing broadband 
%   sounds with own ears) for an exemplary listener (NH12).
%
%      
%
%   See also: mclachlan2021 mclachlan2021 plot_reijniers2014 reijniers2014_metrics
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_mclachlan2021.php


%   #Author : Glen McLachlan (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Settings
%num_dirs = 500; 
%[dirs,~,~,~] = ParticleSampleSphere('N',num_dirs); 
%dirs = load('dirs.mat'); dirs = dirs.dirs;
%[az, el] = cart2sph(dirs(:,1),dirs(:,2),dirs(:,3));

group = 'group_AMTexample2';

%% Get listener's data
SOFA_obj = amt_load('mclachlan2021', 'HRIR_L2354.sofa');

%% Preprocessing source information for both directions
[template, target] = mclachlan2021_preproc(SOFA_obj, group);

%% Run virtual experiments
[doa, params] = mclachlan2021(template, target, group);

%% Calcualte performance measures 
amt_disp('------------------------')
amt_disp('Performance Predictions:')
amt_disp('------------------------')

met = mclachlan2021_metrics(doa);
met.entropy = mean(params.entropy);
met.information = mean(params.information);

%amt_disp(sprintf('Lateral accuracy: %0.2fdeg', met.accL))
amt_disp(sprintf('Lateral RMS: %0.2f', met.rmsL))
%amt_disp(sprintf('Elevation accuracy: %0.2fdeg', met.accE))
amt_disp(sprintf('Polar RMS: %0.2f', met.rmsP))
%amt_disp(sprintf('Quadrant error: %0.6f%%',met.querr))
amt_disp(sprintf('Mean entropy: %0.2fbits', met.entropy))
amt_disp('------------------------')
% 
% %met.nexp = num_exp;
% met.rot_type = rot_type;
% met.rot_size = rot_size;
% met.sig = 1;
% met.coords = target.coords;
% met.post_prob = params.post_prob;

% posterior distribution, this only works if one direction is tested
for i=1:size(params.post_prob,2)
    figure;scatter3(params.template_coords(:,1),params.template_coords(:,2),params.template_coords(:,3),20,squeeze(params.post_prob(:,i,:)),'filled')
    hold on
    plot3(doa.real(1),doa.real(2),doa.real(3),'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12,'LineWidth',2.5);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    colormap(gray)    
    colormap(flipud(gray)) 
end

% this only works if all directions are tested
%plot_reijniers2014(target.coords,params.information);
%title('Information for yaw=0.1 in bits')


