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
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_baumgartner2013.php

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


