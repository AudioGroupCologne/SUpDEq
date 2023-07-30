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
%   2) Data in auxdata/baumgartner2017
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_baumgartner2017.php


%   #Author: Robert Baumgartner (2017)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% TODO: explain Data in description;

%% % Parse input options
definput.flags.type = {'ITE','BTE'};
flags = ltfatarghelper({},definput,varargin);
   
dataPath = 'baumgartner2017';
if flags.do_ITE
  %fn = dir(fullfile(dataPath,'hrtf b_nh*.sofa'));
  ids=[12 14 15 16 17 18 21 22 33 39 41 43 46 53 55 57 58 62 68 71 72];  
  for ii=1:length(ids)
    fn(ii).name=['hrtf b_nh' num2str(ids(ii)) '.sofa'];
  end
elseif flags.do_BTE
  ids=[10];  
  for ii=1:length(ids)
    fn(ii).name=['hrtf b_bte_nh' num2str(ids(ii)) '.sofa'];
  end
end
if isempty(fn)
  error('RB: No files found. Check directories and AMT setup!')
end
data = struct('id',{},'Obj',{});
for ii = 1:length(fn)
  data(ii).id = ids(ii);
  data(ii).Obj = amt_load(dataPath,fn(ii).name);
end

end

