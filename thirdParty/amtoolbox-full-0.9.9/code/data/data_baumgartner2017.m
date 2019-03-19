function data = data_baumgartner2017(varargin)
%DATA_BAUMGARTNER2017  Data from Baumgartner et al. (2017)
%   Usage: data = data_baumgartner2017(flag)
%
%   DATA_BAUMGARTNER2017(flag) returns data from Baumgartner et al. (2017)
%   describing a model for sound externalization.
%
%   The input flag determines the microphone casing and may be
%
%     'ITE'       to obtain a set of 23 standard in-the-ear HRTFs (default)
%     'BTE'       to obtain an exemplary behind-the-ear HRTF
%
%   The fields in the data output contains the following information
%
%     .id         listener ID
%     .Obj        HRTF data in SOFA Format
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2017
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_baumgartner2017.php

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

% TODO: explain Data in description;

%% % Parse input options
definput.flags.type = {'ITE','BTE'};
flags = ltfatarghelper({},definput,varargin);
   
dataPath = fullfile(SOFAdbPath,'baumgartner2017');
if flags.do_ITE
  %fn = dir(fullfile(dataPath,'hrtf b_nh*.sofa'));
  ids=[12 14 15 16 17 18 21 22 33 39 41 43 46 53 55 57 58 62 68 71 72];  
  for ii=1:length(ids)
    fn(ii).name=fullfile(dataPath,['hrtf b_nh' num2str(ids(ii)) '.sofa']);
  end
elseif flags.do_BTE
%   fn = dir(fullfile(dataPath,'hrtf b_bte_nh*.sofa'));
  ids=[10];  
  for ii=1:length(ids)
    fn(ii).name=fullfile(dataPath,['hrtf b_bte_nh' num2str(ids(ii)) '.sofa']);
  end
end
if isempty(fn)
  error('RB: No files found. Check directories and AMT setup!')
end
data = struct('id',{},'Obj',{});
for ii = 1:length(fn)
  data(ii).id = ids(ii);
  data(ii).Obj = SOFAload(fn(ii).name);
end

end
