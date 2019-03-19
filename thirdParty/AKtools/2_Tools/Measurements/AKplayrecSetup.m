% [playDev, recDev, pageSize, fs, maxCHplay, maxCHrec] = AKplayrecSetup(playDev, recDev, pageSize, fs)
%
% GUI based initialization of playrec for audio input/output. Needs playrec
% externals compiled from http://www.playrec.co.uk
%
% See AKmeasureDemo.m for examples
%
% I N P U T / O U T P U T
% all input is optional and only used to choose items in the GUI. If no
% input is specified the default items will be chosen. If you want to skip
% certain input parameters, pass [] or false, e.g. AKplayrecSetup([], 1)
%
% the output returns the parameters selected by the user.
%
% playDev   - playing device id obtained from playrec(getDevices)
% recDev    - recording device id obtained from playrec(getDevices)
% pageSize  - divides audio into blocks of length of pageSize samples
% fs        - sampling rate in Hz
% maxCHplay - number of available channes for the playing device
% maxCHrec  - number of available channes for the recording device
%
% 10/2016  - fabian.brinkmann@tu-berlin.de
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
function [playDev, recDev, pageSize, fs, maxCHplay, maxCHrec] = AKplayrecSetup(playDev, recDev, pageSize, fs)

global Gfs GplayDev GrecDev GpageSize

% --------------------------------------------------- get info from playrec
% get devices
devs.all = playrec('getDevices');

% separate input and output devices
devs.InID    = -1;
devs.InCh    = 0;
devs.InName  = {'No device selected'};

devs.OutID   = devs.InID;
devs.OutCh   = devs.InCh;
devs.OutName = devs.InName;

devs.InCur  = 1;
devs.OutCur = 1;


for nn = 1:numel(devs.all)
    if devs.all(nn).inputChans
        devs.InCur = devs.InCur+1;
        
        devs.InID(devs.InCur, 1)   = devs.all(nn).deviceID;
        devs.InCh(devs.InCur, 1)   = devs.all(nn).inputChans;
        devs.InName{devs.InCur, 1} = ['ID ' num2str(devs.all(nn).deviceID) ': ' devs.all(nn).name ' (' devs.all(nn).hostAPI ', ' num2str(devs.all(nn).inputChans) ' ch.)'];
    end
    if devs.all(nn).outputChans
        devs.OutCur = devs.OutCur+1;
        
        devs.OutID(devs.OutCur, 1)   = devs.all(nn).deviceID;
        devs.OutCh(devs.OutCur, 1)   = devs.all(nn).outputChans;
        devs.OutName{devs.OutCur, 1} = ['ID ' num2str(devs.all(nn).deviceID) ': ' devs.all(nn).name ' (' devs.all(nn).hostAPI ', ' num2str(devs.all(nn).outputChans) ' ch.)'];
    end
end

% -------------------------------------------------- set default parameters 
% available sample rates in Hz
devs.fs = [8000 16000 32000 44100 48000 96000 192000];

if exist('fs', 'var')
    if ~islogical(fs) && ~isempty(fs)
        Gfs = find(devs.fs == fs);
        if isempty(Gfs)
            Gfs = 4;
        end
    else
        Gfs = 4;
    end
else
    Gfs = 4;
end
    
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

if exist('pageSize', 'var')
    if ~islogical(pageSize) && ~isempty(pageSize)
        GpageSize = num2str(pageSize);
    else
        GpageSize = '0';
    end
else
    GpageSize = '0';
end

% --------------------------------------------- display playrec preferences
h.f = AKf(20,13);

a = .05;
b = .12;

set(h.f, 'DockControls', 'off', 'MenuBar', 'none', 'name', 'Audio I/O settings: playrec', 'NumberTitle','off')

uicontrol('Style','text', 'string', 'sampling rate (Hz)', 'units', 'normalized', 'position', [0.12 1-b 0.8 .1], 'fontsize', 16, 'Backgroundcolor', 'w', 'horizontalAlignment', 'left');
h.g1 = uicontrol('Style','popupmenu', 'string', devs.fs, 'units', 'normalized', 'position', [0.1 1-a-b .8 .1], 'fontsize', 16, 'callback', @setFS);
set(h.g1, 'value', Gfs)

uicontrol('Style','text', 'string', 'play device', 'units', 'normalized', 'position', [0.12 1-a-2*b 0.8 .1], 'fontsize', 16, 'Backgroundcolor', 'w', 'horizontalAlignment', 'left');
h.g2 = uicontrol('Style','popupmenu', 'string', devs.OutName, 'units', 'normalized', 'position', [0.1 1-2*a-2*b .8 .1], 'fontsize', 16, 'callback', @setPlayDev);
set(h.g2, 'value', GplayDev)

uicontrol('Style','text', 'string', 'record device', 'units', 'normalized', 'position', [0.12 01-2*a-3*b 0.8 .1], 'fontsize', 16, 'Backgroundcolor', 'w', 'horizontalAlignment', 'left');
h.g3 = uicontrol('Style','popupmenu', 'string', devs.InName, 'units', 'normalized', 'position', [0.1 1-3*a-3*b .8 .1], 'fontsize', 16, 'callback', @setRecDev);
set(h.g3, 'value', GrecDev)

uicontrol('Style','text', 'string', 'page size (audio block length, 0 for default value)', 'units', 'normalized', 'position', [0.12 1-3*a-4*b 0.8 .1], 'fontsize', 16, 'Backgroundcolor', 'w', 'horizontalAlignment', 'left');
h.g4 = uicontrol('Style','edit', 'string', '0', 'units', 'normalized', 'position', [0.1 1-4*a-4*b .8 .1], 'fontsize', 16, 'callback', @setBuffer);
set(h.g4, 'string', GpageSize)

uicontrol('Style','pushbutton', 'string', 'ok', 'units', 'normalized', 'position', [0.625 0.05 .25 .15], 'fontsize', 16, 'callback', 'close gcf');

uiwait(h.f)
pause(.1)

% ------------------------------------------------------------ parse output
maxCHplay  = devs.OutCh(GplayDev);
playDev    = devs.OutID(GplayDev);

maxCHrec   = devs.InCh(GrecDev);
recDev     = devs.InID(GrecDev);

fs         = devs.fs(Gfs);
pageSize = str2double(GpageSize);

% --------------------------------------------- reset playrec if neccessary
if playrec('isInitialised')
    playrec('reset')
end

% -------------------------------------------------------------- initialize
playrec('init', fs, playDev, recDev, maxCHplay, maxCHrec, pageSize)

% print feedback
if playDev < 0 && recDev >=0
    if pageSize
        fprintf('\nplayrec initialized for %s, and %u samples page size\n\n', devs.InName{GrecDev}, pageSize)
    else
        fprintf('\nplayrec initialized for %s, and default page size\n\n', devs.InName{GrecDev})
    end
    warning('AKplayrecSetup:Settings', 'No playback device selected. Audio output not possible!')
elseif recDev < 0 && playDev >= 0
    if pageSize
        fprintf('\nplayrec initialized for %s, and %u samples page size\n\n', devs.InName{GplayDev}, pageSize)
    else
        fprintf('\nplayrec initialized for %s, and default page size\n\n', devs.InName{GplayDev})
    end
    warning('AKplayrecSetup:Settings', 'No record device selected. Audio input not possible!')
elseif playDev < 0 && recDev < 0
    error('AKplayrecSetup:Settings', 'No playback and/or record device selected. Audio In/Out not possible!')
else
    if pageSize
        fprintf('\nplayrec initialized for %s, %s, and %u samples page size\n\n', devs.OutName{GplayDev}, devs.InName{GrecDev}, pageSize)
    else
        fprintf('\nplayrec initialized for %s, %s, and default page size\n\n', devs.OutName{GplayDev}, devs.InName{GrecDev})
    end
end

% clear output and global variables
clear -global Gfs GplayDev GrecDev GpageSize
if ~nargout
    clear playDev recDev pageSize fs maxCHplay maxCHrec
end

end % end of AKplayrecSetup

% ------------------------------------------------------------ ui callbacks
function setFS(varargin)
    global Gfs
    Gfs = get(varargin{1}, 'value');
end

function setPlayDev(varargin)
    global GplayDev
    GplayDev = get(varargin{1}, 'value');
end

function setRecDev(varargin)
    global GrecDev
    GrecDev = get(varargin{1}, 'value');
end

function setBuffer(varargin)
    global GpageSize
    GpageSize = get(varargin{1}, 'string');
end
