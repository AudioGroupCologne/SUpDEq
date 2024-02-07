function varargout=amt_load(model,data,variable)
%AMT_LOAD Load auxiliary data of a model
%   Usage: amt_load(MODEL, DATA);
%
%   AMT_LOAD(model, data) loads the auxiliary data from the file data. The data will loaded
%   from the directory model located in the local auxdata directory given by
%   amt_auxdatapath.
%
%   If the file is not in the local auxdata directory, it will be downloaded from
%   the web address given by amt_auxdataurl.
%
%   For SOFA files, data must be passed in full, i.e. with the '.sofa'
%   extension.
%
%   The following file types are supported:
%
%   - .wav   output will be as that from audioread
%   - .mat   output as that as from load
%   - .sofa  output as that from SOFAload
%   - .csv   output as that from tableread. This is dangerous as it's actual output
%              depends on the Matlab version. Consider storing as mat file.
%   - others   output is the absolute filename
%
%   AMT_LOAD(model, data, variable) loads just a particular variable
%   from the file.
%
%   See also: amt_auxdatapath amt_auxdataurl
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_load.php


%   #Author: Piotr Majdak (2015)
%   #Author: Clara Hollomey (2021)
%   #Author: Clara Hollomey (2023): SOFA support
%   #Author: Piotr Majdak (2023): SOFA support fix, SOFA loaded from auxdata, CSV fixed

% This file is licensed unter the GNU General Public License (GPL) either
% version 3 of the license, or any later version as published by the Free Software
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and
% at <https://www.gnu.org/licenses/gpl-3.0.html>.
% You can redistribute this file and/or modify it under the terms of the GPLv3.
% This file is distributed without any warranty; without even the implied warranty
% of merchantability or fitness for a particular purpose.

[~,kv]=amt_configuration;

localPath=kv.auxdataPath;
localfn=fullfile(localPath,model,data);

% file not found? create directories, and download!
if ~exist(localfn,'file')
    % create dir if not existing
  if ~exist(fullfile(localPath,model),'dir')
    [success,msg]=mkdir(fullfile(localPath,model));
    if success~=1
      error(msg);
    end
  end
    % download
  amt_disp(['Model: ' model '. Downloading auxiliary data: ' data]);
  webfn = [kv.auxdataURL '/' model '/' data];
  webfn(strfind(webfn,'\'))='/';
  webfn=regexprep(webfn,' ','%20');

  if isoctave
    [~, stat]=urlwrite(webfn, localfn);
  else
    try
      outfilename = websave(localfn,webfn);
      stat = 1;
    catch
      stat = 0;
    end
  end

  if ~stat
    for ii = 2:numel(kv.version)
      webfn = num2str([kv.downloadURL kv.version{ii} '/auxdata/' model '/' data]);
      if isoctave
        [~,stat]=urlwrite(webfn,localfn);
      else
        try
          outfilename = websave(localfn,webfn);
          stat = 1;
        catch
          stat = 0;
        end
      end

      if stat
        amt_disp(['Found data in version: ' kv.version{ii}]);
        break;
      end

      if ii == numel(kv.version)
        error(['Unable to download auxdata.']);
      end
    end
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
    if exist('variable','var')
        varargout{1}=load(localfn,variable);
    else
        varargout{1}=load(localfn);
    end
  case '.sofa'
     evalc('output=SOFAload(localfn)');
     varargout{1} = output;
  case '.csv'
     varargout{1}=table2array(readtable(localfn));
  otherwise
    varargout{1}=localfn;
end


