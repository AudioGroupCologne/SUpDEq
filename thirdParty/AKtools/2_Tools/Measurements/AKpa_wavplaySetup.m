% [playDev, recDev, maxCHplay, maxCHrec] = AKplayrecSetup(playDev, recDev)
%
% GUI based initialization of playrec for audio input/output. Needs
% pa_wavplay (included in AKtools)
%
% See AKmeasureDemo.m for examples
%
% I N P U T / O U T P U T
% all input is optional and only used to choose items in the GUI. If no
% input is specified the default items will be chosen. If you want to skip
% certain input parameters, pass [] or false, e.g. AKpa_wavplaySetup([], 1)
%
% the output returns the parameters selected by the user.
%
% playDev   - playing device id obtained from pa_wavplay()
% recDev    - recording device id obtained from pa_wavplay()
% maxCHplay - number of available channes for the playing device
% maxCHrec  - number of available channes for the recording device
%
% 11/2016  - fabian.brinkmann@tu-berlin.de
%

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
function [playDev, recDev, maxCHplay, maxCHrec] = AKpa_wavplaySetup(playDev, recDev)

global GplayDev GrecDev

% ------------------------------------------------ get info from pa_wavplay
% set default device - which is no selection
devs.InID    = -1;
devs.InCh    = 0;
devs.InName  = {'No device selected'};

devs.OutID   = devs.InID;
devs.OutCh   = devs.InCh;
devs.OutName = devs.InName;

devs.InCur  = 1;
devs.OutCur = 1;

% --- get devices ---
try
    devs.all =  evalc('pa_wavplay');
catch ME
    warning('AKpa_wavplaySetup:pa_wavplay', 'pa_wavplay could not be called. Please check if you have ASIO4all installed, and if any soundcard is connected to your computer. For more information see error message below.')
    rethrow(ME);
end

% --- separate input and output devices ---
% find and remove irrelevant output
ID = strfind(devs.all, 'Printing audio devices...');
if ID
    devs.all = devs.all(ID+25:end);
end
devs.all = strtrim(devs.all);
devs.all = regexprep(devs.all, '\n|\r\n|\r', '');

% look for devices
IDall = strfind(devs.all, 'Device');
if ~isempty(IDall)
    
    for nn = 1:numel(IDall)/2
        
        % remove irrelevant output
        devs.all = devs.all(7:end);
        devs.all = strtrim(devs.all);
        
        % get device ID
        ID       = strfind(devs.all, 'Device name:');
        currID   = str2double(devs.all(1:ID(1)-1));
        devs.all = devs.all(ID(1)+12:end);
        devs.all = strtrim(devs.all);
        
        % get device name
        ID       = strfind(devs.all, 'Max Input channels:');
        currName = devs.all(1:ID(1)-1);
        devs.all = devs.all(ID(1)+19:end);
        devs.all = strtrim(devs.all);
        
        % get number of input chanels
        ID       = strfind(devs.all, 'Max Output channels:');
        currIn   = str2double(devs.all(1:ID(1)-2));
        devs.all = devs.all(ID(1)+20:end);
        devs.all = strtrim(devs.all);
        
        %get number of output channels
        ID = strfind(devs.all, 'Device');
        if ID
            currOut  = str2double(devs.all(1:ID(1)-1));
            % BUGFIX: old line produced an error on WIN10 with Scarlett 18i20 1st gen. sound card
            % Old line: devs.all = devs.all(ID(1)-1:end);
            devs.all = devs.all(ID(1):end);
            devs.all = strtrim(devs.all);
        else
            currOut = str2double(devs.all(1:end));
        end
        
        % write to struct
        if currIn
            devs.InCur = devs.InCur+1;
            
            devs.InID(devs.InCur, 1)   = currID;
            devs.InCh(devs.InCur, 1)   = currIn;
            devs.InName{devs.InCur, 1} = ['ID ' num2str(currID) ': ' currName ' (' num2str(currIn) ' ch.)'];
        end
        if currOut
            devs.OutCur = devs.OutCur+1;
            
            devs.OutID(devs.OutCur, 1)   = currID;
            devs.OutCh(devs.OutCur, 1)   = currOut;
            devs.OutName{devs.OutCur, 1} = ['ID ' num2str(currID) ': ' currName ' (' num2str(currOut) ' ch.)'];
        end
    end
end

clear ID IDall currID currName currIn currOut nn

% -------------------------------------------------- set default parameters 
if exist('playDev', 'var')
    if ~islogical(playDev) && ~isempty(playDev)
        GplayDev = find(devs.OutID == playDev);
    else
        GplayDev = 1;
    end
else
    GplayDev = 1;
end

if exist('recDev', 'var')
    if ~islogical(recDev) && ~isempty(recDev)
        GrecDev = find(devs.InID == recDev);
    else
        GrecDev = 1;
    end
else
    GrecDev = 1;
end

% --------------------------------------------- display playrec preferences
h.f = AKf(20,13);

a = .05;
b = .12;

set(h.f, 'DockControls', 'off', 'MenuBar', 'none', 'name', 'Audio I/O settings: pa_wavplay', 'NumberTitle','off')

uicontrol('Style','text', 'string', 'play device', 'units', 'normalized', 'position', [0.12 1-a-2*b 0.8 .1], 'fontsize', 16, 'Backgroundcolor', 'w', 'horizontalAlignment', 'left');
h.g2 = uicontrol('Style','popupmenu', 'string', devs.OutName, 'units', 'normalized', 'position', [0.1 1-2*a-2*b .8 .1], 'fontsize', 16, 'callback', @setPlayDev);
set(h.g2, 'value', GplayDev)

uicontrol('Style','text', 'string', 'record device', 'units', 'normalized', 'position', [0.12 01-2*a-3*b 0.8 .1], 'fontsize', 16, 'Backgroundcolor', 'w', 'horizontalAlignment', 'left');
h.g3 = uicontrol('Style','popupmenu', 'string', devs.InName, 'units', 'normalized', 'position', [0.1 1-3*a-3*b .8 .1], 'fontsize', 16, 'callback', @setRecDev);
set(h.g3, 'value', GrecDev)

uicontrol('Style','pushbutton', 'string', 'ok', 'units', 'normalized', 'position', [0.625 0.05 .25 .15], 'fontsize', 16, 'callback', 'close gcf');

uiwait(h.f)
pause(.1)

% ------------------------------------------------------------ parse output
maxCHplay  = devs.OutCh(GplayDev);
playDev    = devs.OutID(GplayDev);

maxCHrec   = devs.InCh(GrecDev);
recDev     = devs.InID(GrecDev);

% clear output and global variables
clear -global GplayDev GrecDev
if ~nargout
    clear playDev recDev maxCHplay maxCHrec
end

end % end of AKplayrecSetup

% ------------------------------------------------------------ ui callbacks
function setPlayDev(varargin)
    global GplayDev
    GplayDev = get(varargin{1}, 'value');
end

function setRecDev(varargin)
    global GrecDev
    GrecDev = get(varargin{1}, 'value');
end
