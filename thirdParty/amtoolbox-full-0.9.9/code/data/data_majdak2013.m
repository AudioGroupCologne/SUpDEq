function data = data_majdak2013(varargin)
%DATA_MAJDAK2013 Listener specific localization in saggital planes
%   Usage: data = data_majdak2013(condition)
%
%   Output parameters:
%     data.id    : listener ID
%     data.mtx   : experimental data matrix conaining 9 colums
%                  col 1: target azimuth
%                  col 2: target elevation
%                  col 3: response azimuth
%                  col 4: response elevation
%                  col 5: lateral angle of target
%                  col 6: polar angle of target
%                  col 7: lateral angle of response
%                  col 8: polar angle of response
%
%   DATA_MAJDAK2013(condition) returns listener-specific experimental data
%   from Majdak et al.  (2013) testing localization performance in sagittal
%   planes for low-pass filtered and spectrally warped DTFs.
%
%   The condition flag may be one of:
%
%     'BB'   Broadband DTFs (baseline condition). This is the default.
%     'LP'   Low-pass filtered (at 8.5kHz) DTFs
%     'W'    Spectrally warped (2.8-16kHz warped to 2.8-8.5kHz) DTFs
%
%   References:
%     P. Majdak, T. Walder, and B. Laback. Effect of long-term training on
%     sound localization performance with spectrally warped and band-limited
%     head-related transfer functions. J. Acoust. Soc. Am., 134:2148-2159,
%     2013.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_majdak2013.php

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

% AUTHOR: Robert Baumgartner

%% Check input options

% Define input flags
definput.flags.condition = {'BB','LP','W'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Extract data
x=amt_load('majdak2013','data.mat');

C = find(ismember(x.condition,flags.condition));

for ll = 1:length(x.subject)
  
  data(ll).mtx = x.subject(ll).expData{C}(:,1:8);
  data(ll).id = x.subject(ll).id;

end
