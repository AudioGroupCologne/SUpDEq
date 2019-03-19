function op1=amt_version(varargin)
%amt_version Help on the AMToolbox
%   Usage:  amt_version;
%           v=amt_version('version');
%           mlist=amt_version('modules');
%
%   AMT_VERSION displays some general help on the AMT.
%
%   AMT_VERSION('version') returns the version number.
%
%   AMT_VERSION('modules') returns a cell array of installed modules and
%   corresponding version numbers.
%
%   See also:  amt_start
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/amt_version.php

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

%   AUTHOR : Peter Søndergaard. 
%   MODIFICATIONS: 09.07.2017 Piotr Majdak

bp=amt_basepath;

definput.keyvals.versiondata=[];
definput.keyvals.modulesdata=[];
definput.flags.mode={'general','version','modules','authors'};

[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_general
  amt_disp(' ');
  amt_disp('--- AMT - The Auditory Modeling Toolbox. ---');
  amt_disp(' ')

  amt_disp(['Version ',kv.versiondata]);
  amt_disp(' ');
  amt_disp('Installed modules:');
  amt_disp(' ');
  amt_disp('Name:            Version:  Description');
  modinfo=amt_version('modules');
  for ii=1:length(modinfo);
    s=sprintf(' %-15s %7s  %s',modinfo{ii}.name,modinfo{ii}.version, ...
	      modinfo{ii}.description);
    amt_disp(s);
  end;

end;
  
if flags.do_version
  op1=kv.versiondata;
end;

if flags.do_modules
  op1={};
  for ii=1:numel(kv.modulesdata)
    
    p=kv.modulesdata{ii};
    
    % Get the first line of the help file
    [FID, MSG] = fopen ([bp,p.name,filesep,'Contents.m'],'r');
    if FID==-1
      error('Module %s does not contain a Contents.m file.',p.name);
    end;
    firstline = fgetl (FID);
    fclose(FID);
    
    
    % Load the information into the cell array.	
    op1{ii}.name=p.name;
    op1{ii}.version=p.version;
    op1{ii}.description=firstline(2:end);
  end;
end;

