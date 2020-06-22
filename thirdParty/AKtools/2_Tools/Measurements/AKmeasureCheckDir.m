% This script is called by AKmeasureDemo from AKtools
% See AKmeasureDemo.m for examples

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

% check if base directory exists
function data = AKmeasureCheckDir(data)

global makeDirs

if isempty(data.dir)
    data.dir  = uigetdir(fullfile(fileparts(which('AKtoolsStart'))), 'Location for Saving Measurements');
end

if data.dir ~= 0
    
    %ask to create the direcory if it does not exist
    if ~exist(data.dir, 'dir')
        
        h.f = AKf(20,7);
        set(h.f, 'DockControls', 'off', 'MenuBar', 'none', 'name', 'Enter file name (files will be saved to data.dir)', 'NumberTitle','off')
        h.e = uicontrol('Style','text', 'string', {data.dir  ' does not exist. Do you want to create it?'}, 'units', 'normalized', 'position', [0.05 0.4 0.9 .5], 'fontsize', 16, 'Backgroundcolor', 'w');
        uicontrol('Style','pushbutton', 'string', 'yes', 'units', 'normalized', 'position', [0.25 0.1 .25 .2], 'fontsize', 16, 'callback', {@createDirs, data});
        uicontrol('Style','pushbutton', 'string', 'no', 'units', 'normalized', 'position', [.5 0.1 .25 .2], 'fontsize', 16, 'callback', @createDirs);
        
        uiwait(h.f)
        pause(.1)
    else
        makeDirs = true;
        disp('Output directory is valid.')
    end
    
    % make (sub)directories
    if makeDirs
        
        if ~exist(data.dir, 'dir')
            mkdir(data.dir)
        end
        
        if ~exist(fullfile(data.dir, 'Plots'), 'dir')
            mkdir(fullfile(data.dir, 'Plots'))
        end
        
        if ~exist(fullfile(data.dir, 'Data'), 'dir')
            mkdir(fullfile(data.dir, 'Data'))
        end
        
        if ~exist(fullfile(data.dir, 'Reference'), 'dir')
            mkdir(fullfile(data.dir, 'Reference'))
        end
    end
else
    warning('AKmeasure:SaveDir', 'Please select a directory for saving the data')
end

clearvars -global makeDirs

end

function createDirs(varargin)
global makeDirs

if strcmpi(get(varargin{1}, 'string'), 'yes')
    makeDirs = true;
    disp('Output directory created.')
else
    makeDirs = false;
    disp('No output directory created - try again.')
end

close gcf

end
