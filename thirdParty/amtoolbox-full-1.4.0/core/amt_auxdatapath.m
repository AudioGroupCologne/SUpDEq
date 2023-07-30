function auxPath=amt_auxdatapath(newPath)
%AMT_AUXDATAPATH Local path to the auxiliary data
%   Usage: auxpath=amt_auxdatapath
%          amt_auxdatapath(newpath)
%
%   auxPath=AMT_AUXDATAPATH returns the path of the directory containing
%   auxiliary data.
%
%   Default path to the auxiliary data is amt_basepath/auxdata.
% 
%   AMT_AUXDATAPATH(newpath) sets the path of the directory for further calls
%   of AMT_AUXDATAPATH.
%
%   See also: amt_auxdataurl amt_load amt_basepath
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_auxdatapath.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   #Author: Piotr Majdak (2015)
%   #Author: Clara Hollomey (2021)

mlock;
persistent CachedPath;
caller = dbstack;
last = numel(caller);

if exist('newPath','var')
  CachedPath=newPath;
elseif isempty(CachedPath)
  CachedPath=fullfile(amt_basepath, 'auxdata');
end
auxPath=CachedPath;

if ~strcmp('amt_configuration', caller(last).name)
  amt_configuration('auxdataPath', auxPath);
end
  


