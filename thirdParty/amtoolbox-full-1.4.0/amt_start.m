function amt_start(varargin)
%amt_start   Start the Auditory Modeling Toolbox (AMT)
%   Usage:  amt_start;
%           amt_start(flags);
%
%   AMT_START starts the AMT. This command must be
%   run before using any of the function in the AMT.
%
%   AMT_START('install') queries the user for yes/no to download and install
%   all available third-party toolboxes. Then, it executes amt_mex to compile 
%   the binaries on the system. 
%
%
%   Cache:
%   ------
%
%   AMT uses cache to store precalculated results because some of the AMT functions
%   require a large processing time. Depending on the machine and the model, it might take
%   even days. The global cache mode is controlled on start-up of the AMT. To change the
%   global cache mode choose a flags:
%
%     'normal'      Use cached package as far as possible. This is default.
%                   This is kind of demonstration mode and very convenient
%                   for fast access of results like plotting figures.
%                   This option, however, may by-pass the actual processing and thus
%                   does not always test the actual functionality of a model.
%                   If the cached package locally not available, downloaded from the internet.
%                   If remotely not available, enforce recalculation.
%
%     'cached'      Enforce to use cached package. If the cached package is
%                   locally not available, it will be downloaded from the internet.
%                   If it is remotely not available, an error will be thrown.
%
%     'redo'        Enforce the recalculation of the package. This option
%                   actually tests the calculations.
%
%     'localonly'   Package will be recalculated when locally
%                   not available. Do not connect to the internet.
%
%   Many AMT functions support the cache mode as input flag in order to
%   overwrite the global cache mode. See AMT_CACHE for more details.
%
%
%   Auxiliary data
%   --------------
%
%   Most of the models require auxiliary data. The AMT will download these data on-demand.
%   The download URL for the auxiliary data is given by amt_auxdataurl.
%   The target directory for the auxiliary data is given by amt_auxdatapath.
%   If you want to run the AMT offline, download the auxiliary data first.
%
%   Some of the auxiliary data are HRTFs. The AMT will download the HRTFs on-demand.
%   The download URL for the HRTFs is given by SOFAdbURL.
%   The target directory for the HRTFs is given by SOFAdbPath.
%   If you want to run the AMT offline, download the HRTFs first.
%
%   Output
%   ------
%
%   The output of the messages to the command line can be controlled by one
%   of the following flags:
%
%     'verbose'        All output will be displayed. This is default.
%
%     'documentation'  starts the AMT in the documentation compiling
%                      mode. The output of calculation progress will be suppressed.
%
%     'silent'         All output will be suppressed.
%
%
%   Go to http://www.amtoolbox.org/doc.php for the full documentation.
%
%   See also:  amt_mex amt_load amt_cache amt_disp
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/amt_start.php


%   #Author : Peter L. Soendergaard (2013) 
%   #Author: Piotr Majdak (2013-)
%   #Author: Clara Hollomey(2020-2023)
%   #Author: Piotr Majdak (2023): fix for using cache flags in documentation mode

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Start AMT====================================================================

%where are we?
resetpath = pwd;
basepath = which('amt_start.m');
basepath = basepath(1:end-numel('amt_start.m'));
%store the path as it was before the addition of basepath
%(ltfatstart deletes all global variables, so it can only
%be stored as global afterwards)
P = addpath(basepath);

%how should I display?
if any(strcmp(varargin, 'silent'))
    dispFlag = 'silent';
    silent = 1;
    documentation = 0;
elseif any(strcmp(varargin, 'documentation'))
    dispFlag = 'documentation';
    silent = 0;
    documentation = 1;
else
    dispFlag = 'verbose';
    silent = 0;
    documentation = 0;
end

if any(strcmp(varargin, 'install'))
    install = 1;
    debug = 1;
    dispFlag = 'verbose';
else
    install = 0;
    debug = 0;
end


% if isempty(varargin) || any(strcmp(varargin, 'install'))
%     dispFlag = 'verbose';
% end
% 
% if any(strcmp(varargin, 'verbose')) || isempty(varargin)
%     debug = 0;
%     dispFlag = 'verbose';
% else
%     debug = 0;
% end
% 
% if any(strcmp(varargin, 'silent'))
%     silent = 1;
%     dispFlag = 'silent';
%     debug = 0;
% else
%     silent = 0;
%     debug = 0;
% end
% 
% if any(strcmp(varargin, 'documentation'))
%     documentation = 1;
%     dispFlag = 'documentation';
%     debug = 0;
% else
%     documentation = 0;
%     debug = 0;
% end
% 
% if any(strcmp(varargin, 'install'))
%     install = 1;
%     debug = 1;
%     dispFlag = 'verbose';
% else
%     install = 0;
%     debug = 0;
% end

if ~exist('arg_amt_configuration','file')
   cd(fullfile(basepath,'defaults'));
end

%now, that the core folder is available, load the default configuration...
%(to set the configuration, we need ltfatarghelper, which we do not yet have)
definput = arg_amt_configuration;

if exist('OCTAVE_VERSION','builtin')
  definput.keyvals.interpreter = 'Octave';
  definput.keyvals.interpreterversion = version;
  warning('off','Octave:shadowed-function'); % remove warnings about functions shadowing other functions
else
  definput.keyvals.interpreter = 'Matlab';
  definput.keyvals.interpreterversion = version;  
end

if isunix
  definput.keyvals.system = 'Linux';
elseif ismac
  definput.keyvals.system = 'Mac';
elseif ispc
  definput.keyvals.system = 'Windows';
else
  definput.keyvals.system = 'other';
end

if strcmp(definput.keyvals.interpreter, 'Octave')
  expectedFolders = numel(definput.keyvals.amt_folders);
else
  expectedFolders = numel(definput.keyvals.amt_folders) - 1;
end

for ii = 1:expectedFolders
    definput.keyvals.amt_folders(ii).name = [basepath, definput.keyvals.amt_folders(ii).name];
    addpath(definput.keyvals.amt_folders(ii).name);
end

if ~silent
  disp(' ');
  disp('****************************************************************************');
  disp(' ');
  disp(['The Auditory Modeling Toolbox (AMT) version ', definput.keyvals.version{1}(5:end),'.']);
  disp('Brought to you by the AMT Team. See http://amtoolbox.org for more information.');
  disp(' ');
  disp('****************************************************************************');
  disp(' ');
end
cd(basepath);
definput.keyvals.amtrunning = 1;
%% Starting the LTFAT

%amt_subdir is now available...so let's search for LTFAT
thirdpartypath=fullfile(basepath,'thirdparty');
addpath(thirdpartypath);

ltfatpath = amt_installtoolboxes('ltfat', 'ltfatstart.m', definput.keyvals.ltfat, thirdpartypath, 1);  
addpath(ltfatpath);

if silent, x='ltfatstart(''nojava'',0);'; else x='ltfatstart(''nojava'');'; end
S = evalc(x);
if debug
  disp(S); 
else
  if ~silent, disp(['* LTFAT: ' S(1:min([find(S==10) length(S)+1])-1)]); end
end

% Check for the correct LTFAT version.
s=ltfathelp('version');
s_r='2.0.0'; % set the required version
v=sscanf(s,'%d.%d.%d'); v(4)=0;
v_r=sscanf(s_r,'%d.%d.%d');

if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || ...
   (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) )
    
   error(['You need LTFAT >= ' s_r ' to work with AMT. ' ...
             'Update your package from ' definput.keyvals.ltfat '.']);
end
definput.keyvals.ltfatfound = ltfatpath;
definput.keyvals.ltfatrunning = 1;

%% check if all necessary toolboxes (MATLAB) and packages (Octave) are installed.

if isoctave
  [signaltb, loaded] = amt_checkoctpkg('signal');
  if signaltb
    if ~silent 
      if loaded
      disp('* Signal package: loaded.');
      definput.keyvals.signal = 1;
      elseif ~loaded
        disp('* Signal package: not found. Load this package and then re-run amt_start.');
        definput.keyvals.signal = 0;
      end
    end
  else
    if ~silent 
      disp('* Signal package: not found. Load this package and then re-run amt_start.');
      definput.keyvals.signal = 0;
    end
  end
  
  [statisticstb, loaded] = amt_checkoctpkg('statistics');
  if statisticstb
    if ~silent 
      if loaded
      disp('* Statistics package: loaded.');
      definput.keyvals.statistics = 1;
     elseif ~loaded  
      disp('* Statistics package: not found. Load this package and then re-run amt_start.');
      definput.keyvals.statistics = 0;
     end
   end
  else
    if ~silent 
      disp('* Statistics package: not found. Load this package and re-run amt_start.'); 
      definput.keyvals.statistics = 0;
    end
  end  
  
  [netcdftb, loaded] = amt_checkoctpkg('netcdf');
  if netcdftb
    if ~silent
      if loaded
        disp('* Netcdf package: loaded.');
        definput.keyvals.netcdf = 1;
      elseif ~loaded
        disp('* Netcdf package: not found. Load this package and then re-run amt_start.');
        definput.keyvals.netcdf = 0;
      end
    end
  else
    if ~silent 
      disp('* Netcdf package: not found. Load this package and re-run amt_start.');
      definput.keyvals.netcdf = 0;
    end
  end
  [optimtb, loaded] = amt_checkoctpkg('optim');
  if optimtb
    if ~silent
      if loaded
        disp('* Optim package: loaded.');
        definput.keyvals.optim = 1;
      elseif ~loaded
        disp('* Optim package: not found. Load this package and then re-run amt_start.');
        definput.keyvals.optim = 0;
      end
    end
  else
    if ~silent 
      disp('* Optim package: not found. Load this package and re-run amt_start.');
      definput.keyvals.optim = 0;
    end
  end   
else
  installedToolboxes = ver;
  definput.keyvals.optim = 'NA';
  definput.keyvals.netcdf = 'NA';
  for ii = 1:numel(installedToolboxes)
      if strcmp(installedToolboxes(ii).Name, 'Signal Processing Toolbox')
          if ~silent, disp('* Signal Processing Toolbox: Found.'); end
          definput.keyvals.signal = 1;
      elseif strcmp(installedToolboxes(ii).Name, 'Statistics and Machine Learning Toolbox')
          if ~silent, disp('* Statistics and Machine Learning Toolbox: Found.'); end
          definput.keyvals.statistics = 1;
      end
  end
    
  if ~isfield(definput.keyvals, 'signal')
      if ~silent, disp('* Signal Processing Toolbox: Not found.'); end
  elseif ~isfield(definput.keyvals, 'statistics')
      if ~silent, disp('* Statistics and Machine Learning Toolbox: Not found.'); end
  end
    
end


%% now, check for the optional thirdparty toolboxes
if ~silent,
  disp(' ');
  disp('****************************************************************************');
  disp(' ');
  disp('Starting optional third-party toolboxes:');
end
%% SOFA package=================================================================
% Search for SOFA package
targetfile = 'SOFAstart.m';
sofapath = [];
if exist(targetfile,'file')
  %check if the targetfile is already in the path
  sofapath = fileparts(which(targetfile)); 
else
  %check if the targetfile exists in the installationpath
  sofapath = amt_subdir(fullfile(thirdpartypath, targetfile));
  if ~isempty(sofapath),
    addpath(sofapath(1).name(1:end-numel(targetfile)));
  end
end

if install && isempty(sofapath)
  sofapath = amt_installtoolboxes('sofa', 'SOFAstart.m', definput.keyvals.sofa, thirdpartypath, 0);
  if ~isempty(sofapath)
    addpath(sofapath);
  end
end 

if exist('SOFAstart','file')
 
  if silent, x='SOFAstart(''restart'');'; else x='SOFAstart(''restart'');'; end
  S = evalc(x);
  if debug, 
    disp(S); 
  else
    if ~silent, disp(['* SOFA: ' S(1:min([find(S==10) length(S)+1])-1)]); end
  end;
  
  SOFAdbPath(fullfile(basepath,'hrtf'));
  %SOFAdbURL('https://www.sofacoustics.org/data');
  SOFAdbURL(definput.keyvals.hrtfURL);
  % SOFAdbURL is a default path and will be overwritten later
  warning('off','SOFA:upgrade'); % disable warning on upgrading older SOFA files
  warning('off','SOFA:load'); % disable warnings on loading SOFA files
  definput.keyvals.sofarunning = 1;
  definput.keyvals.sofafound = sofapath;
else
  definput.keyvals.sofarunning = 0;
  definput.keyvals.sofafound = 'NO';
  if ~silent, disp('* SOFA: Toolbox not available. No SOFA support, limited HRTF support only. Run amt_start(''install''); to change that.');  end
end


%% SFS package==================================================================
targetfile = 'SFS_start.m';
sfspath = [];
if exist(targetfile,'file')
  %check if the targetfile is already in the path
  sfspath = fileparts(which(targetfile));
  addpath(sfspath);
else
  %check if the targetfile exists in the installationpath
  sfspath = amt_subdir(fullfile(thirdpartypath, targetfile));
  if ~isempty(sfspath),
    addpath(sfspath(1).name(1:end-numel(targetfile)));
  end
end

if install && isempty(sfspath)
  sfspath = amt_installtoolboxes('sfs', 'SFS_start.m', definput.keyvals.sfs, thirdpartypath, 0);
  if ~isempty(sfspath)
    addpath(sfspath);
  end
end

if exist('SFS_start','file')
%if we have sfs installed, we should delete rms.m because of syntax conflicts
sfsdeletepath=fileparts(which('SFS_start.m'));
if isoctave
  if exist(fullfile(sfsdeletepath,'SFS_octave','rms.m'),'file')
	  delete(fullfile(sfsdeletepath,'SFS_octave','rms.m'));
  end
else
  if exist(fullfile(sfsdeletepath,'SFS_general','rms.m'),'file')
	  delete(fullfile(sfsdeletepath,'SFS_general','rms.m'));
  end
end
end
% Start

if exist('SFS_start','file')
  
  S = evalc('SFS_start');
  s=SFS_version;
  if ~isempty(S),      
    if ~silent, disp(['* SFS: ' S(1:min([find(S==10) length(S)+1])-1)]); end
  else
    if ~silent, 
      disp(['* SFS: Toolbox version ' s ' loaded']); 
        if isoctave, warning('off','Octave:num-to-str'); end %Octave only: turn off warning about sfs info display 
    end    
  end;
    % check required version
  s_r='2.5.0'; % Required version
  v=sscanf(s,'%d.%d.%d'); v(4)=0;
  v_r=sscanf(s_r,'%d.%d.%d');    
  if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
      error(['You need SFS >= ' s_r ' to work with AMT. ' ...
        'Please update your package from https://github.com/sfstoolbox/sfs ']);
  end    
  
  definput.keyvals.sfsrunning = 1;
  definput.keyvals.sfsfound = sfspath;
else
  definput.keyvals.sfsrunning = 0;
  definput.keyvals.sfsfound = 'NO';
  if ~silent, disp('* SFS: SFS Toolbox NOT available. Run amt_start(''install''); to change that.');  end
end

%% Circular Statistics Toolbox (Directional Statistics)=========================
targetfile = 'kuipertable.mat';
circstatpath = [];
if exist(targetfile,'file')
  %check if the targetfile is already in the path
  circstatpath = fileparts(which(targetfile));
  if ~isempty(circstatpath)
    addpath(circstatpath);
  end
else
  %check if the targetfile exists in the installationpath
  circstatpath = amt_subdir(fullfile(thirdpartypath, targetfile));
  if ~isempty(circstatpath) 
      addpath(circstatpath(1).name(1:end-numel(targetfile)));
  end
end

if install && isempty(circstatpath)
  circstatpath = amt_installtoolboxes('circstat', 'kuipertable.mat', definput.keyvals.circstat, thirdpartypath, 0);
  if ~isempty(circstatpath)
    addpath(circstatpath);
  end
end

if ~isempty(circstatpath)
  definput.keyvals.circstatrunning = 1;
  definput.keyvals.circstatfound = circstatpath;
  if ~silent, disp('* CIRCSTAT: Circular statistics toolbox found.'); end
else
  definput.keyvals.circstatrunning = 0;
  definput.keyvals.circstatfound = 'NO';
  if ~silent, disp('* CIRCSTAT: Circular statistics toolbox NOT available. Run amt_start(''install''); to change that.'); end
end
cd(basepath);

%% binauralSH Toolbox======================================================
targetfile = 'binauralSH_start.m';
binauralshpath = [];
if exist(targetfile,'file')
  %check if the targetfile is already in the path
  binauralshpath = fileparts(which(targetfile));
  addpath(binauralshpath);
else
  %check if the targetfile exists in the installationpath
  binauralshpath = amt_subdir(fullfile(thirdpartypath, targetfile));
  if ~isempty(binauralshpath)
    addpath(binauralshpath(1).name(1:end-numel(targetfile)));
  end
end

if install && isempty(binauralshpath)
  binauralshpath = amt_installtoolboxes('binauralSH', 'binauralSH_start.m', definput.keyvals.binauralSH, thirdpartypath, 0);
  if ~isempty(binauralshpath)
    addpath(binauralshpath);
  end
end


if exist('binauralSH_start','file')
  S = evalc('binauralSH_start');
  if debug, 
    disp(S); 
  elseif isempty(S),
    if ~silent, disp('* BINAURALSH: BinauralSH toolbox loaded.'); end
  else
    if ~silent, disp(['* BINAURALSH: ' S(1:min([find(S==10) length(S)+1])-1)]); end
  end;  
  definput.keyvals.binshrunning = 1;
  definput.keyvals.binshfound = binauralshpath;
  if ~silent
      if isoctave
        warning('off','Octave:num-to-str');%Octave only: turn off warning about sfs info display
      end 
  end
else
  definput.keyvals.binshrunning = 0;
  definput.keyvals.binshfound = 'NO';
  if ~silent, disp('* BINAURALSH: BinauralSH toolbox NOT available. Run amt_start(''install''); to change that.'); end
end
cd(resetpath);

if ~silent
  disp(' ');
  disp('****************************************************************************');
  disp(' ');
  disp('Internal configuration:');
end

%% Set the correct path to remote HRTFs
if exist('SOFAdbURL','file'),
     SOFAdbURL(definput.keyvals.hrtfURL);
end
%% Initialize aux data, cache, and display starting information
if ~silent, 
  disp('  ');
  disp(['  Auxiliary data (local): ' definput.keyvals.auxdataPath]);
  disp(['  Auxiliary data (web): ' definput.keyvals.auxdataURL]);
end

%check python version
[status ,result]=(system('python --version'));
definput.keyvals.python = result;
[~,kv]=amt_configuration(definput);

%set the cache and disp for later
cacheFlag = 'normal';
if any(strcmp(varargin, 'redo'))
    cacheFlag = 'redo';
elseif any(strcmp(varargin, 'localonly'))
    cacheFlag = 'localonly';
elseif any(strcmp(varargin, 'cached'))
    cacheFlag = 'cached';
elseif documentation
    cacheFlag = 'redo';
end

[flags,~]=amt_configuration(dispFlag, cacheFlag);
global userdefaultpath
userdefaultpath = P;

if ~silent, 
  switch flags.cachemode
    case 'normal'
      disp('  Cache mode: Download precalculated results. Examples:');
      disp('              exp_model(...)        shows precalculated results');
      disp('              exp_model(...,''redo'') enforces recalculation');
    case 'localonly'
      disp('  Cache mode: Use local cache or recalculate. Do not connect to remote cache.');
    case 'cached'
      disp('  Cache mode: Use cache or throw error. Do not recalcalculate.');
    case 'redo'
      disp('  Cache mode: Recalculate always. Be patient!');
  end
  disp(' ');
  disp('  Check your configuration via [flags, keyvalues] = amt_configuration().');
  disp('Type "help amt_start" for more details...');
end
if install
    try
    amt_mex;
    catch
        amt_disp('mex files were not compiled.');
    end
end

function filepath = amt_installtoolboxes(toolboxname, targetfile, downloadpath, installationpath, required)
%this is a helper function to install thirdparty (matlab) toolboxes to the AMT.

if exist(targetfile,'file')
  %check if the targetfile is already in the path
  filepath = fileparts(which(targetfile));  
else
  %check if the targetfile exists in the installationpath
  filepath = amt_subdir(fullfile(installationpath, targetfile));

  if ~isempty(filepath)
    %if found, great, let's add it to the path
    filepath = filepath.name;
    filepath = filepath(1:end-(numel(targetfile)+1));
  else
    %else, ask the user if it should be installed 
    disp(' ');
    disp([upper(toolboxname) ' package is neither in the current path nor was it found in the thirdparty folder.']);
    if required
      disp(['Unable to continue without ' upper(toolboxname) '.']);
    end
    if ~required
      toolboxInstall ='n';
      %while toolboxInstall == 'n'
        toolboxInstall = input('Download it to ./thirdparty? (y/n)>','s');
      %end
    else
      toolboxInstall = 'y';
    end
    
    if strcmp(toolboxInstall, 'y')

        disp('Downloading...');
        %download from the path in the configuration     
        downloadTarget = installationpath;
        webfn = downloadpath;
        archiveName = sprintf('%sZIP',toolboxname);
        tokenfn = fullfile(downloadTarget, archiveName);
        if ~exist(downloadTarget, 'dir') mkdir(downloadTarget); end
  
        if exist('OCTAVE_VERSION','builtin')
          [~, stat]=urlwrite(webfn, tokenfn);
        else
          if verLessThan('matlab','9.2')           
              options = weboptions('Timeout',1000);
              outfilename = websave(tokenfn,webfn);
          else
              o = weboptions();
              certificateFile = o.CertificateFilename; % was here before 1.2.1-dev
              %options = weboptions('Timeout',100, 'CertificateFilename','');
              options = weboptions('CertificateFilename','');  % was here before 1.2.1-dev
              options=weboptions; options.CertificateFilename=(''); % PM for 1.2.1-dev
              outfilename = websave(tokenfn,webfn,options); % PM for 1.2.1-dev
          end
%           outfilename = websave(tokenfn,webfn);
          if exist(outfilename, 'file')
              stat = 1;
          end
        end
        
        if ~stat
          warning(['Unable to download file from remote: ' webfn]);
          disp('Please check your internet settings.');
        else
          disp('Unzipping...');
          downloaded=unzip(fullfile(installationpath, archiveName), fullfile(installationpath,toolboxname));
          if exist(fullfile(installationpath, [archiveName, '.zip']))
            delete(sprintf('%s',fullfile(installationpath, [archiveName, '.zip'])));
          elseif exist(fullfile(installationpath, archiveName))
            delete(sprintf('%s',fullfile(installationpath, archiveName)));
          end
        end

        %now, check again if the targetfile can be found  
        filepath = amt_subdir(fullfile(installationpath, toolboxname, targetfile));   
      
        if ~isempty(filepath)
        %if the targetfile can be found now...
          filepath = filepath.name;
          filepath = filepath(1:end-(numel(targetfile)+1));        
          disp([upper(toolboxname) ' installation was successful, it now resides within ./thirdparty.']);
          disp(' ');
        else
          filepath = [];
          if required
            error([toolboxname ' installation was unsuccessful. If you did not get a download warning, ' 10 ...
                 'a zip file may still reside within ./thirdparty - please unzip it and manually run amt_start again. ' 10 ...
                 'If not, please check your internet settings and/or copy the toolbox manually to ./thirdparty.' 10 ...
                 'Then run amt_start again.']);
          else 
            warning([toolboxname ' installation was unsuccessful. If you did not get a download warning, ' 10 ...
                 'a zip file may still reside within ./thirdparty - please unzip it and manually run amt_start again. ' 10 ...
                 'If not, please check your internet settings and/or copy the toolbox manually to ./thirdparty. ' 10 ...
                 'Then run amt_start again.']);
          end      
        end
      
    else %this is just to totally ensure that there either was an error or that the download succeeded 
        filepath = []; 
        if required 
          error([toolboxname ' not found. Unable to continue without it.' 10 ...
             '   Download the toolbox from' downloadpath ' and' 10 ...
             '   copy it to amtoolbox/thirdparty/. Then run amt_start again.']);
        end
           
             
    end
  end
end


function [installed, loaded] = amt_checkoctpkg(pkgname)
%this function checks for the installed and loaded Octave packages
installed = 0;
loaded = 0;

[~,info]=pkg('list');
for ii = 1:numel(info)
  if strcmp(pkgname, info{ii}.name)
    installed = 1;
    if ~info{ii}.loaded
      try
        eval(sprintf('pkg load %s', pkgname));
        loaded = 1;
      catch
        loaded = 0;
      end
    else
      loaded = 1;
    end
  end    
end    
    


