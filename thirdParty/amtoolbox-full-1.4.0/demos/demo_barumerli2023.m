%DEMO_BARUMERLI2023 Demo for full sphere localization model from Barumerli et al. (2023)
%
%   DEMO_BARUMERLI2023(flag) demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   on the full sphere using the localization model from Barumerli et al. (2023).
%
%   Figure 1: Baseline prediction
% 
%      This demo computes the baseline prediction (localizing broadband 
%      sounds with own ears) for an exemplary listener (NH12).
%
%      Averaged polar and lateral accuracy
%
%   See also: barumerli2023 exp_barumerli2023 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_barumerli2023.php


%   #Author: Roberto Barumerli (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


amt_disp('Loading real data from Majdak et al. 2010')

data_majdak = data_majdak2010('Learn_M');
data_majdak = data_majdak(6);

%% real subject
% print metrics computed on real data 
real = barumerli2023_metrics(data_majdak.mtx, 'middle_metrics');

amt_disp('Metrics from actual data')
amt_disp(sprintf('lateral_bias [deg]:\t%0.2f', real.accL))
amt_disp(sprintf('lateral_rms_error [deg]:\t%0.2f', real.rmsL))
amt_disp(sprintf('elevation_bias [deg]:\t%0.2f', real.accP))
amt_disp(sprintf('local_rms_polar [deg]:\t%0.2f', real.rmsP))
amt_disp(sprintf('quadrant_err [%%]:\t%0.2f', real.querr))

%% simulated subject
% load HRTF dataset
amt_disp('Model simulation starting...')

sofa_obj = amt_load('barumerli2023','ARI_NH12_hrtf_M_dtf 256.sofa');

% compute features for target and templates 
amt_disp(sprintf('\nGenerating feature space over the uniformily sampled sphere'))
[template, target] = barumerli2023_featureextraction(sofa_obj, 'pge');

% estimate each direction in the target struct twice
amt_disp('Simulating subject answers')
m = barumerli2023(...
                'template', template, ...
                'target', target, ...
                'num_exp', 2, ...
                'sigma_itd', 0.569, ...
                'sigma_ild', 0.7, ...
                'sigma_spectral', 1.25, ...
                'sigma_prior', 11.5, ...
                'sigma_motor', 11);

% compute metrics over estimated directions
sim = barumerli2023_metrics(m, 'middle_metrics');

amt_disp('Metrics from simulated data')
amt_disp(sprintf('lateral_bias [deg]:\t%0.2f', sim.accL))
amt_disp(sprintf('lateral_rms_error [deg]:\t%0.2f', sim.rmsL))
amt_disp(sprintf('elevation_bias [deg]:\t%0.2f', sim.accP))
amt_disp(sprintf('local_rms_polar [deg]:\t%0.2f', sim.rmsP))
amt_disp(sprintf('quadrant_err [%%]:\t%0.2f', sim.querr))

%% ESTIMATE ONE DIRECTION
amt_disp(sprintf('\nSimulating single estimation'))

targ_az = 0;
targ_el = 45;

% compute features only for the single bianural sound
% (the template has been computed above)
target_single = barumerli2023_featureextraction(sofa_obj, 'pge', 'target','targ_az', targ_az, 'targ_el', targ_el);

% simulate single estimation process
[m, doa, doa_real] = barumerli2023(...
                'template', template, ...
                'target', target_single, ...
                'num_exp', 1, ...
                'sigma_itd', 0.569, ...
                'sigma_ild', 0.7, ...
                'sigma_spectral', 1.25, ...
                'sigma_prior', 11.5, ...
                'sigma_motor', 11);

% print resulted estimation
for i=1:size(targ_az,1)
    amt_disp(sprintf('target [az,el]: [%.2f, %.2f]', targ_az(i), targ_el(i)))
    amt_disp(sprintf('response [az,el]: [%.2f, %.2f]\n', m(i, 3), m(i, 4)))
end

plot_reijniers2014(template.coords.return_positions('cartesian'), ...
                    doa.posterior, ...
                    'target', doa_real.return_positions('cartesian'));
title('Posterior distribution and target direction')

