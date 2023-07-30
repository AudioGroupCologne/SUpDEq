function data = data_majdak2010(varargin)
%DATA_MAJDAK2010 Listener specific localization performance
%   Usage: data = data_majdak2010(condition)
%
%   DATA_MAJDAK2010(condition) returns listener-specific experimental data
%   from Majdak et al. (2010) testing localization performance for various
%   experimental methods.
%
%   The data struct comprises the following fields:
%
%     'id'    listener ID
%
%     'mtx'   experimental data matrix containing 9 columns
%
%             - target azimuth
%             - target elevation
%             - response azimuth
%             - response elevation
%             - lateral angle of target
%             - polar angle of target
%             - lateral angle of response
%             - polar angle of response
%
% 
%   The condition flag may be one of:
%
%     'HMD_M'   Head-mounted display and manual pointing. Testing of naive
%               subjects.
%     'HMD_H'   Head-mounted display, head pointing, naive subjects.
%     'Dark_M'  Dark room, manual pointing, naive subjects.
%     'Dark_H'  Dark room, head pointing, naive subjects.
%     'Learn_M' Acoustic learning condition with manual pointing. This is the default.
%     'Learn_H' Acoustic learning condition with head pointing.
%
%   References:
%     P. Majdak, M. J. Goupell, and B. Laback. 3-D localization of virtual
%     sound sources: Effects of visual environment, pointing method and
%     training. Atten Percept Psycho, 72:454--469, 2010.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_majdak2010.php


%   #Author: Robert Baumgartner (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Check input options

% Define input flags
definput.flags.condition = {'Learn_M','Learn_H','HMD_M','HMD_H','Dark_M','Dark_H'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Extract data
x=amt_load('majdak2010','data.mat');

C = find(ismember(x.condition,flags.condition));

for ll = 1:length(x.subject)
  
  if not(isempty(x.subject(ll).expData{C}))
    data(ll).mtx = real(x.subject(ll).expData{C}(:,1:8));
  end
  data(ll).id = x.subject(ll).id;

end


