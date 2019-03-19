%DEMO_BAUMGARTNER2014 Demo for sagittal-plane localization model from Baumgartner et al. (2014)
%
%   DEMO_BAUMGARTNER2014(flag) demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   for a listener of the listener pool and the median plane using the 
%   sagittal-plane localization model from Baumgartner et al. (2014).
%
%   Figure 1: Baseline prediction
% 
%      This demo computes the baseline prediction (localizing broadband 
%      sounds with own ears) for an exemplary listener (NH58).
%
%      Predicted polar response angle probability of subject NH58 as a  
%      function of the polar target angle with probabilities encoded by
%      brigthness.
%
%   See also: baumgartner2014 exp_baumgartner2014 baumgartner2014_virtualexp
%   localizationerror
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_baumgartner2014.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR : Robert Baumgartner

%% Settings

subID = 'NH58';   % subject ID of exemplary listener
lat = 0;          % lateral target angle in degrees
runs = 3;         % # of virtual experimental runs
   

%% Get listener's data

s = data_baumgartner2014('pool');   % load data of listener pool
ids = find(ismember({s.id},subID));  % index of exemplary listener


%% Run model with individual sensitivity S

[p,rang,tang] = baumgartner2014(s(ids).Obj,s(ids).Obj,'S',s(ids).S,'lat',lat);


%% Run virtual experiment

m = baumgartner2014_virtualexp(p,tang,rang,'runs',2);


%% Calcualte performance measures 

amt_disp('Performance Predictions:')
amt_disp('------------------------')

% via expectancy values:
[qe,pe] = baumgartner2014_pmv2ppp(p,tang,rang,'print'); 

% and/or via responses drawn from virtual experiments
[f,r] = localizationerror(m,'sirpMacpherson2000');
perMacpherson2003 = localizationerror(m,f,r,'perMacpherson2003');
amt_disp(['Local polar error rate (%)        ' num2str(perMacpherson2003,'%4.1f')])


%% Plot results

figure;
plot_baumgartner2014(p,tang,rang,m(:,6),m(:,8));
title(['Baseline prediction for ' s(ids).id]);

