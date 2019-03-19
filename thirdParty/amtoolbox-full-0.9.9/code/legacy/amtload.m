function varargout=amtload(model,data,variable)

warning('Warning: AMTLOAD will be removed in a future release. Use AMT_LOAD instead. ');  

% based on amt_load from 9.7.2017
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/amtload.php

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

localfn=fullfile(amt_auxdatapath,model,data);
  % file not found? create directories, and download!
if ~exist(localfn,'file')
    % create dir if not existing
  if ~exist(fullfile(amt_auxdatapath,model),'dir'), 
    [success,msg]=mkdir(fullfile(amt_auxdatapath,model));
    if success~=1
      error(msg);
    end
  end
    % download
  amt_disp(['Model: ' model '. Downloading auxiliary data: ' data],'progress');
  webfn=[amt_auxdataurl '/' model '/' data];
  webfn(strfind(webfn,'\'))='/';
  webfn=regexprep(webfn,' ','%20');        
  [~,stat]=urlwrite(webfn,localfn);
  if ~stat
    error(['Unable to download file: ' webfn]);
  end
end
  % load the content
[~,~,ext] = fileparts(localfn);
switch lower(ext)
  case '.wav'
    [y,fs]=audioread(localfn);
    varargout{1}=y;
    varargout{2}=fs;
  case '.mat'
    if exist('variable','var'),
        varargout{1}=load(localfn,variable);
    else
        varargout{1}=load(localfn);
    end
  otherwise
    varargout{1}=localfn;
end
