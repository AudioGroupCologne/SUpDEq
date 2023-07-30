%DEMO_BAUMGARTNER2016 Demo for sagittal-plane localization model from Baumgartner et al. (2016)
%
%   DEMO_BAUMGARTNER2016 demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   for a listener of the listener pool and the median plane using the 
%   sagittal-plane localization model from Baumgartner et al. (2016).
%
%   Figure 1: Baseline prediction
% 
%      This demo computes the baseline prediction (localizing broadband 
%      sounds with own ears) for an exemplary listener with outer-hair-cell dysfunction.
%
%   See also: baumgartner2016 exp_baumgartner2016 baumgartner2014_virtualexp
%   localizationerror
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_baumgartner2016.php


%   #Author: Robert Baumgartner (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Settings

condition = 'baseline';
subID = 'NH33';   % subject ID of exemplary listener
lat = 0;          % lateral target angle in degrees
SPL = 60;         % SPL of target sound
runs = 1;         % # of virtual experimental runs
cohc = 0.4;       % outer-hair-cell dysfunction
do_plot = true;   % flag for plotting predicted response PMVs
performanceMetric = 'querrMiddlebrooks';
   
%% Get listener's data

s = data_baumgartner2016;   % load data of listener pool
ids = find(ismember({s.id},subID));  % index of exemplary listener

%% Run model

[err,pred,m] = baumgartner2016(s(ids).Obj,s(ids).Obj,'ID',subID,'Condition',condition,... 
  'lat',lat,'S',s(ids).S,'SPL',SPL,'cohc',cohc,performanceMetric);

%% Comparison with actual performance metric 

err_exp = localizationerror(s(ids).itemlist,performanceMetric);

disp(['Actual normal-hearing performance (',performanceMetric,'): ',num2str(err_exp,'%2.1f')])
disp(['Predicted performance (',performanceMetric,'): ',num2str(err,'%2.1f')])

%% Plot results

figure;
plot_baumgartner2014(pred.p,pred.tang,pred.rang);
title([subID,', C_{OHC} = ',num2str(cohc,'%1.1f')])

