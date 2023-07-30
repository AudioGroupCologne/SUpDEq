function definput = arg_amt_configuration(definput)


%parameters to be set at runtime (not to be changed by the user)================
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_amt_configuration.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
bp = which('amt_start');

definput.keyvals.path = num2str(bp(1:end-numel('amt_start.m')-1));

definput.keyvals.system = 'not determined';
definput.keyvals.interpreter = 'not determined';
definput.keyvals.interpreterversion = 'not determined';
definput.keyvals.python = 'not determined';


%thirdparty toolboxes
definput.keyvals.ltfatfound = 'not initialized';%these will be paths after initialization...
definput.keyvals.sofafound = 'not initialized';
definput.keyvals.sfsfound = 'not initialized';
definput.keyvals.circstatfound = 'not initialized';
definput.keyvals.binshfound = 'not initialized';

definput.keyvals.ltfatrunning = 0;
definput.keyvals.sofarunning = 0;
definput.keyvals.sfsrunning = 0;
definput.keyvals.circstatrunning = 0;
definput.keyvals.binshrunning = 0;
definput.keyvals.signal = 0;  
definput.keyvals.statistics = 0;
definput.keyvals.optim = 0;  
definput.keyvals.netcdf = 0;  
definput.keyvals.amtrunning = 0;  

%amt folders to install=========================================================
definput.keyvals.amt_folders(1).name = 'defaults';%this folder always needs to be the first in the list
definput.keyvals.amt_folders(2).name = 'core';
definput.keyvals.amt_folders(3).name = 'data';
definput.keyvals.amt_folders(4).name = 'common';
definput.keyvals.amt_folders(5).name = 'experiments';
definput.keyvals.amt_folders(6).name = 'environments';
definput.keyvals.amt_folders(7).name = 'legacy';
definput.keyvals.amt_folders(8).name = 'mex';
definput.keyvals.amt_folders(9).name = 'models';
definput.keyvals.amt_folders(10).name = 'modelstages';
definput.keyvals.amt_folders(11).name = 'plot';
definput.keyvals.amt_folders(12).name = 'signals';
definput.keyvals.amt_folders(13).name = 'thirdparty';
definput.keyvals.amt_folders(14).name = 'demos';
definput.keyvals.amt_folders(15).name = 'oct';%this folder always needs to be the last in the list


%version-dependent parameters (not to be changed)===============================
%definput.keyvals.core = num2str(fullfile(definput.keyvals.path, 'general'));
definput.keyvals.default = num2str(fullfile(definput.keyvals.path, 'defaults'));
definput.keyvals.ltfat = num2str('https://github.com/ltfat/ltfat/archive/refs/tags/v2.5.0.zip');
definput.keyvals.sofa = num2str('https://sourceforge.net/projects/sofacoustics/files/SOFA%20Toolbox%20for%20Matlab%20and%20Octave%202.1.4.zip/download');
definput.keyvals.sfs = num2str('https://github.com/sfstoolbox/sfs-matlab/archive/2.5.0.zip');
definput.keyvals.circstat = num2str('https://de.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/10676/versions/22/download/zip');
definput.keyvals.binauralSH = num2str('https://github.com/isaacengel/BinauralSH/archive/refs/heads/main.zip');

definput.keyvals.downloadURL = 'http://amtoolbox.org/';
curver=['amt-' cell2mat(textread(fullfile(definput.keyvals.path,'amt_version'),'%s'))];
definput.keyvals.version = {curver, 'amt-1.3.0', 'amt-1.2.0', 'amt-1.1.0', 'amt-1.0.0'};
definput.keyvals.hrtfURL = num2str([definput.keyvals.downloadURL definput.keyvals.version{1} '/hrtf']);

%parameters to be changed by the user===========================================
definput.keyvals.cacheURL = num2str([definput.keyvals.downloadURL definput.keyvals.version{1} '/cache']);
definput.keyvals.auxdataURL = num2str([definput.keyvals.downloadURL definput.keyvals.version{1} '/auxdata']);
definput.keyvals.auxdataPath = num2str([definput.keyvals.path, filesep, 'auxdata']);
definput.flags.cachemode={'global','normal','cached','localonly','redo'};
definput.flags.disp={'verbose','documentation','silent'};

