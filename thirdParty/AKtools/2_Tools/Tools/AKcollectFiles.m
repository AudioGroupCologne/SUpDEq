% [fileList, nameList] = AKcollectFiles(baseDir, fileName)
% searches all subfolders of baseDir to find files specified by fileName
%
% For example
% AKcollectFiles('/Users/A_Custic/audio, '*.wav')
% searches for wav files in /audio and all it's subfolders
%
% I N P U T
% baseDir  - the complete path to the folder that should be searched
% fileName - expression to search for, e.g. '*.wav' will return all wave
%            files and '*raw*.wav' will return all wave files whose name
%            contains the word 'raw'. By default, all files are returned,
%            which is equivalent to passing '*'. Multiple expressions can
%            be passed in a struct, e.g., {'*.wav', '*.mat'}
%
% O U T P U T
% fileList - a cell array that lists all found files by full path and file
%            name
% nameList - a cell array that lists only the file names without the full
%            path
%
%
% fabian.brinkmann@tu-berlin.de
% 02/2018   - initial code

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [fileList, nameList] = AKcollectFiles(baseDir, fileName)

% check if multiple fileNames are passed
if nargin == 2 && iscell(fileName)
    fileList = [];
    nameList = [];
    
    for ff = 1:numel(fileName)
        [fileListTmp, nameListTmp] = AKcollectFiles(baseDir, fileName{ff});
        
        fileList = [fileList; fileListTmp]; %#ok<AGROW>
        nameList = [nameList; nameListTmp]; %#ok<AGROW>
    end
    return
end

% allocate space
fileList  = cell(1e6, 1);
nameList  = fileList;

% default input
if nargin == 1
    fileName = '*';
end

dirs = regexp(genpath(baseDir),'[^:;]*','match')';

i = 1;

for nn = 1:numel(dirs)
    % get all filenames in current directory
    curFiles = dir( fullfile( dirs{nn}, fileName ) );
    curFiles = struct2cell(curFiles)';
    id       = cell2mat(curFiles(:,5));
    curFiles = curFiles(~id,1:2);
    
    % discard files that start with . and ~
    id = ~startsWith(curFiles(:,1), '.') & ~startsWith(curFiles(:,1), '~');
    curFiles = curFiles(id,:);
    
    % get full path
    curNames = curFiles(:,1);
    curFiles = strcat( curFiles(:,2), filesep, curFiles(:,1) );
    % add to output
    fileList(i:i+numel(curFiles)-1) = curFiles;
    nameList(i:i+numel(curFiles)-1) = curNames;
    
    % increase counter
    i = i + numel(curFiles);
end

fileList = fileList(1:i-1);
nameList = nameList(1:i-1);