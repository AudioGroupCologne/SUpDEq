%DEMO_BAUMGARTNER2013 Demo for sagittal-plane localization model from Baumgartner et al. (2013)
%
%   DEMO_BAUMGARTNER2013(flag) demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   for a listener of the listener pool and the median SP using the 
%   sagittal-plane localization model from Baumgartner et al. (2013).
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
%   See also: baumgartner2013
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_baumgartner2013.php


%   #Author: Robert Baumgartner (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Settings

flag='NH58';  % subject ID
lateral = 0;  % lateral target angle in degrees
   

%% Search subject

s = data_baumgartner2013('pool');
for ids = 1:length(s)
    if strcmp(flag,s(ids).id)
        break
    end
end


%% Run model

[targets,tang] = extractsp(lateral,s(ids).Obj);
[p,rang] = baumgartner2013(targets,s(ids).Obj,'u',s(ids).u,'lat',lateral);

figure;
plot_baumgartner2013(p,tang,rang);
title(['Baseline prediction for ' s(ids).id]);

[qe,pe,pb] = baumgartner2013_pmv2ppp(p,tang,rang,'print');



