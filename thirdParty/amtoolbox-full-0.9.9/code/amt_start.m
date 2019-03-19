function amt_start(varargin)
%amt_start   Start the Auditory Modeling Toolbox
%   Usage:  amt_start;
%           amt_start(flags);
%
%   AMT_START starts the AMT. This command must be
%   run before using any of the function in the AMT.
%
%   Requirements to run the AMT
%   ---------------------------
%
%   1) Linear Time Frequency Analysis Toolbox (*LTFAT*): Download LTFAT from
%      <http://ltfat.sourceforge.net/> and unpack to the prepared
%      directory thirdparty/ltfat. Alternatively, save the LTFAT anywhere and add
%      the root LTFAT (just the root!) path to the search path.
%
%   2) The SOFA API (version >= 1.0): Download the SOFA API from <http://sourceforge.net/projects/sofacoustics>
%      and unpack to the prepared directory thirdparty/SOFA. Alternatively, save
%      the SOFA API anywhere and add the root path of the SOFA API (just the root!) to the search path.
%
%   3) The SFS Toolbox (version >= 2.4.0). Download from <https://github.com/sfstoolbox/sfs-matlab/releases> and
%      unpack to the prepared directory thirdparty/sfs. Alternatively, save
%      the SFS Toolbox anywhere and add the path to the search path.
%
%   4) Python (version >= 2.6 but < 3) with the packages numpy and scipi.
%      On Linux, type sudo apt-get install python-scipy python-numpy.
%      On Windows, intall python from <https://www.python.org/>, add python.exe to the Windows search path,
%      and install the packages separately.
%   5) AMT_MEX executed without any errors. See AMT_MEX for compiler requirements.
%
%   Note that some models may further require other toolboxes. See the corresponding documentation for more details.
%
%
%   Cache:
%   ------
%
%   AMT uses cache to store precalculated results because some of the AMT functions
%   require large processing time. Depending on the machine and the model, it might take
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
%   Go to http://amtoolbox.sourceforge.net/doc.php for the full documentation.
%
%   See also:  amt_mex amt_load amt_cache
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/amt_start.php

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


%   AUTHOR : Peter L. Soendergaard (2013), Piotr Majdak (2013-)


%% Start AMT
bp=amt_basepath;

% Search for LTAFT package
if ~exist('ltfatstart','file')
  ltfatpath=fullfile(bp,'thirdparty','ltfat');
  if exist(ltfatpath,'dir')
    addpath(ltfatpath);
  end
end

% Load the version number
[FID, MSG] = fopen ([bp,'amt_version'],'r');
if FID == -1
    error(MSG);
else
    version = fgetl (FID);
    fclose(FID);
end

% Check if 'silent' present in the flags
silent=0;
if exist('OCTAVE_VERSION','builtin'), args=argv; else args=varargin; end
 for ii=1:numel(args)
   s=lower(args{ii});
   if strcmp(s,'silent') || strcmp(s,'-q')
     silent=1;
   end;
 end;


if ~silent
  disp('  ');
  disp(['AMT version ',version,'. (C) Piotr Majdak and the AMT team.']);
  disp(['See http://amtoolbox.sourceforge.net for more information.']);
  disp('  ');
  disp('Starting toolboxes...');
end;

%% LTFAT package

% Start LTFAT
% if ~silent, disp('*** Starting LTFAT ***'); end
if exist('ltfatstart','file')
  curdir=pwd; % save the current directory
  if silent, ltfatstart('nojava',0); else ltfatstart('nojava'); end;
  cd(curdir); % go back to the saved directory because LTFAT 2.1.1 is messing up the directory in Octave
else
  error(['LTFAT package could not be found. Unable to continue.' 10 ...
        'Download LTFAT from http://ltfat.sourceforge.net ' 10 ...
        'and copy to amtoolbox/thirdparty/ltfat.']);
end

% Check for the correct version.
s=ltfathelp('version');
s_r='2.0.0'; % set the required version
v=sscanf(s,'%d.%d.%d'); v(4)=0;
v_r=sscanf(s_r,'%d.%d.%d');
if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
    error(['You need LTFAT >= ' s_r ' to work with AMT. ' ...
      'Please update your package from http://ltfat.sourceforge.net ']);
end

%% SOFA package

% Search for SOFA package
basepath=which('amt_start');
basepath=basepath(1:end-11);
if ~exist('SOFAstart','file')
  sofapath=fullfile(basepath,'thirdparty','SOFA','API_MO');
  if exist(sofapath,'dir')
    addpath(sofapath);
  end
end

% Start SOFA
if exist('SOFAstart','file')
  SOFAdbPath(fullfile(basepath,'hrtf'));
  SOFAdbURL('http://www.sofacoustics.org/data'); % This is a default path and will be overwritten later
  if silent, SOFAstart('silent'); else SOFAstart('short'); end
	warning('off','SOFA:upgrade');	% disable warning on upgrading older SOFA files
	warning('off','SOFA:load'); % disable warnings on loading SOFA files
else
  if ~silent,
  disp('SOFA package could not be found. Continue without SOFA support.');
  disp(['For SOFA support please download the package ' ...
        'from http://sofacoustics.sourceforge.net ' ...
        'and copy to amtoolbox/thirdparty/SOFA.']);
  end
end

%% SFS package

% Search for the package
basepath=which('amt_start');
basepath=basepath(1:end-11);
if ~exist('SFS_start','file')
  sfspath=fullfile(basepath,'thirdparty','sfs');
  if exist(sfspath,'dir')
    addpath(sfspath);
  end
end

% Delete rms.m from the SFS package because of naming conflict
sfspath=fileparts(which('SFS_start.m'));
if exist(fullfile(sfspath,'SFS_general','rms.m'),'file'),
	delete(fullfile(sfspath,'SFS_general','rms.m'));
end

% Start
% if ~silent, disp('*** Starting SFS ***'); end
if exist('SFS_start','file')
  SFS_start;
  s=SFS_version; s_r='2.4.0'; % set the required version
  if ~silent
      disp(['Sound Field Synthesis Toolbox, version ' s ...
            '. For help, see http://matlab.sfstoolbox.org.']);
  end
  v=sscanf(s,'%d.%d.%d'); v(4)=0;
  v_r=sscanf(s_r,'%d.%d.%d');
  if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
      error(['You need SFS >= ' s_r ' to work with AMT. ' ...
        'Please update your package from https://github.com/sfstoolbox/sfs ']);
  end

elseif ~silent,
  disp('SFS package could not be found. Continue without SFS support.');
  disp(['For SFS support please download the package ' ...
        'from https://github.com/sfstoolbox/sfs ' ...
        'and copy to amtoolbox/thirdparty/sfs.']);
end

%% Install AMT modules
% A directory called DIRNAME containing a file 'DIRNAMEinit.m' is
% considered as a module.
% DIRNAMEinit.m must set the variable 'status' with the following value:
%  0: disabled module, don't add to the search path
%  >0: add to the search path.
%  2: add th amt_version as a module.

% add root of the AMT to the path
addpath(basepath);

modules={};
nplug=0;

% List all files in base directory
d=dir(basepath);

for ii=1:length(d)
  if d(ii).isdir
    if ~(d(ii).name(1)=='.')
      % The file is a directory and it does not start with '.' This could be a module
      name=d(ii).name;
      if exist([bp,name,filesep,name,'init.m'],'file')
          % Set 'status' to zero if the module forgets to define it.
        status=0;
        module_version=version;
        addpath([bp,name]);

        eval([name,'init']);
        if status>0
          if status==1
            nplug=nplug+1;
            modules{nplug}.name=name;
            modules{nplug}.version=module_version;
          end;
        else
          rmpath([bp,name]);
        end;
      end;

    end;
  end;
end;

%% define default start-up behaviour
flags=amt_flags(varargin); % amt_disp and other amt-related functions work now!

% As comp is now in the path, we can set defaults for ltfatarghelper
ltfatsetdefaults('amt_version','versiondata',version,...
                 'modulesdata',modules);

%% Set the correct path to remote HRTFs
if exist('SOFAdbURL','file'),
    SOFAdbURL(['http://www.sofacoustics.org/data/amt-' amt_version('version') '/hrtf']);
end

%% Initialize aux data, cache, and display starting information
amt_disp('  ');
amt_disp('AMT configuration:');
amt_disp(['  Auxiliary data (local): ' amt_auxdatapath]);
amt_disp(['  Auxiliary data (web): ' amt_auxdataurl]);
if strcmp(flags.cachemode,'global'), flags.cachemode='normal'; end
amt_cache('setMode',flags.cachemode);
switch flags.cachemode
  case 'normal'
    amt_disp('  Cache mode: Download precalculated results. Examples:');
    amt_disp('              exp_model(...)        shows precalculated results');
    amt_disp('              exp_model(...,''redo'') enforces recalculation');
  case 'localonly'
    amt_disp('  Cache mode: Use local cache or recalculate. Do not connect to remote cache.');
  case 'cached'
    amt_disp('  Cache mode: Use cache or throw error. Do not recalcalculate.');
  case 'redo'
    amt_disp('  Cache mode: Recalculate always (be patient!).');
end
amt_disp(' ');
amt_disp('Type "help amt_start" for more details...');


