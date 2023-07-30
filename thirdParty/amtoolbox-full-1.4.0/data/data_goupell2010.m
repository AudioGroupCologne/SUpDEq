function data = data_goupell2010(varargin)
%DATA_GOUPELL2010 Localization performance in sagittal planes
%   Usage: data = data_goupell2010(condition)
%          data = data_goupell2010(lat, dlat, condition)
%
%   Output parameters:
%     data : structure
%
%   The condition flag may be one of:
%     'BB'   Broadband DTFs (baseline condition). This is the default.
%     'CL'   Click trains with unlimited number of channels
%     'N24'  24 vocoder channels
%     'N18'  18 vocoder channels
%     'N12'  12 vocoder channels
%     'N9'   9 vocoder channels
%     'N6'   6 vocoder channels
%     'N3'   3 vocoder channels
%
%
%   The 'data' struct contains the following fields:
%     '.id'    : listener ID
%     '.mtx'   : experimental data matrix conaining 9 colums
%
%
%   The columns contain:
%     'col1'  target azimuth
%     'col2'  target elevation
%     'col3'  response azimuth
%     'col4'  response elevation
%     'col5'  lateral angle of target
%     'col6'  polar angle of target
%     'col7'  lateral angle of response
%     'col8'  polar angle of response
%
%   Listener-specific experimental data from Goupell et al. (2010) testing
%   localization performance in sagittal planes for various numbers of
%   channels of a GET vocoder.
%
%   References:
%     M. J. Goupell, P. Majdak, and B. Laback. Median-plane sound
%     localization as a function of the number of spectral channels using a
%     channel vocoder. J. Acoust. Soc. Am., 127:990--1001, 2010.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_goupell2010.php


%   #Author: Robert Baumgartner

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Check input options

% Define input flags
definput.flags.condition = {'BB','CL','N24','N18','N12','N9','N6','N3'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Extract data
x=amt_load('goupell2010','data.mat');

C = find(ismember(x.condition,flags.condition));

for ll = 1:length(x.subject)
  
  data(ll).mtx = x.subject(ll).expData{C}(:,1:8);
  data(ll).id = x.subject(ll).id;

end


