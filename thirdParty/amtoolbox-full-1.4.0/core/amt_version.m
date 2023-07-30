function op1=amt_version(varargin)
%AMT_VERSION Help on the AMToolbox
%   Usage:  amt_version;
%           v=amt_version('version');
%
%   AMT_VERSION displays some general AMT banner.
%
%   AMT_VERSION('version') returns the version number.
%
%
%   See also:  amt_start
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_version.php


%   #Author : Peter Søndergaard. 
%   #Author : Piotr Majdak (2017)  
%   #Author : Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

bp=amt_basepath;
[flags, kv] = amt_configuration;
definput.keyvals.versiondata=kv.version{1};
definput.keyvals.modulesdata=[];
definput.flags.mode={'general','version','modules','authors'};

[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_general
  amt_disp(' ');
  amt_disp('--- AMT - The Auditory Modeling Toolbox. ---');
  amt_disp(' ')

  amt_disp(['-------------Version ',kv.versiondata,'-------------']);
  amt_disp(' ');


end;
  
if flags.do_version
  op1=kv.versiondata;
end;

if flags.do_modules
  amt_disp('This functionality has been deprecated.');
end;



