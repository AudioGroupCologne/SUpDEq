function data = data_majdak2013ctc(varargin)
%DATA_MAJDAK2013CTC Listener-specific localization in sagittal planes
%   Usage: data = data_majdak2013ctc(condition)
%          data = data_majdak2013ctc(lat, dlat, condition)
%
%   DATA_MAJDAK2013CTC(condition) returns listener-specific experimental
%   data from Majdak et al. (2013) testing localization performance in
%   sagittal planes for repeated HRTF measurements motivated by CTC binaural
%   synthesis.
%
%   The data struct comprises the following fields:
%
%     'id'     listener ID
%     'mtx'    experimental data matrix containing 9 columns
%              - target azimuth
%              - target elevation
%              - response azimuth
%              - response elevation
%              - lateral angle of target
%              - polar angle of target
%              - lateral angle of response
%              - polar angle of response
%
%
%   The condition flag may be one of:
%
%     'Learn' Last 300 trials of acoustical training with visual feedback.
%     'A'     First HRTF measurement. This is the default.
%     'B'     Second HRTF measurement.
%
%
%   References:
%     P. Majdak, B. Masiero, and J. Fels. Sound localization in
%     individualized and non-individualized crosstalk cancellation systems.
%     J. Acoust. Soc. Am., 133:2055--68, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_majdak2013ctc.php


%   #Author: Robert Baumgartner (2015)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Check input options

% Define input flags
definput.flags.condition = {'A','B','Learn'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Extract data
x=amt_load('majdak2013ctc','data.mat');

C = find(ismember(x.condition,flags.condition));

for ll = 1:length(x.subject)
  
  data(ll).mtx = real(x.subject(ll).expData{C}(:,1:8));
  data(ll).id = x.subject(ll).id;

end


