function auxURL=amt_auxdataurl(newURL)
%AMT_AUXDATAURL URL of the auxiliary data
%   Usage: url=amt_auxdataurl
%          amt_auxdataurl(newurl)
%
%   url=AMT_AUXDATAURL returns the URL of the web address containing
%   auxiliary data.
% 
%   AMT_AUXDATAURL(newurl) sets the URL of the web address for further calls
%   of AMT_AUXDATAURL.
%
%   See also: amt_auxdatapath amt_load
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_auxdataurl.php


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
persistent AuxDataURL;
caller = dbstack;
last = numel(caller);

[~, kv] = amt_configuration;

if exist('newURL','var')
  AuxDataURL=newURL;
elseif isempty(AuxDataURL)
  AuxDataURL = kv.auxdataURL;
end
auxURL=AuxDataURL;

if ~strcmp('amt_configuration', caller(last).name)
  amt_configuration('auxdataURL', auxURL);
end


