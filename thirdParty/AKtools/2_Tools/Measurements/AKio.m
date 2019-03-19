% [recording, playDevice, recordDevice, engine] = AKio(x, playCH, recordCH, fs, playDevice, recordDevice, engine, pageSize)
% can be used to play and record audio with AKtools. For this purpose
% external audio engines (pa_wavplay, and playrec) are used. Playrec can be
% used from Windows and Mac, hoever might produce glitches in Windows.
% pa_wavplay can be used from Windows, only.
%
% For audio playback pass an audio signal
% AKio(x),
% for recording specify the number of samples to record
% AKio(44100);
% and for playback and recording specify input and output channel(s)
% AKio(x, 1, 1);
%
% see AKioDemo.m for more examples
% see AKmeasureDemo.m for a more sophisticated use case
%
% I N P U T:
% x            - audio signal of size [N C] if you want to playback audio
%                or integer N specifiying the record duration in samples.
%                (N = number of samples; C = number of channels). If x is
%                an audio signal the record duration is set to size(x,1)
% playCH       - vector specifying the playback channels e.g. [1 3], or 2
%                (default is 1:size(x,2))
% recordCH     - vector specifying the record channels e.g. [1 3], or 2
%                (default is 1)
% fs           - sampling rate in Hz (default = 44100)
% playDevice   - ID of playDevice obtained from AKpa_wavplaySetup.m or
%                AKplayrecSetup.m. If this is not passed or empty, the
%                playDevice can be chosen from a GUI
% recordDevice - ID of recordDevice obtained from AKpa_wavplaySetup.m or
%                AKplayrecSetup.m. If this is not passed or empty, the
%                recordDevice can be chosen from a GUI
% engine       - 'playrec' (default for Mac), or
%                'pa_wavplay' (default for Windows)
% pageSize     - if engine = 'playrec' this specifiey the number of samples
%                per channel that are passed to playrec in each block.
%                (playrec takes care of this by default, but you can
%                experiment with it in case of glitches)
%
% O U T P U T:
% recording    - recorded audio signal of size [N numel(recordCH)] where
%                the first is recorded from recordCH(1) etc.
% playDevice   - the used playDevice (see input parameters)
% recordDevice - the used recordDevice (see input parameters)
% engine       - the used audio engine (see input parameters)
%
% 12/2016 - fabian.brinkmann@tu-berlin.de

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
% limitations under  the License.
function [recording, playDevice, recordDevice, engine] = AKio(x, playCH, recordCH, fs, playDevice, recordDevice, engine, pageSize)

% ---------------------------------------------------------- 1. parse input
if numel(x) == 1
    x = zeros(x,1);
end
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('playCH', 'var')
    playCH = [];
end
if ~exist('recordCH', 'var')
    recordCH = [];
end
if ~exist('playDevice', 'var')
    playDevice = [];
end
if ~exist('recordDevice', 'var')
    recordDevice = [];
end
if ~exist('engine', 'var')
    if isunix
        engine = 'playrec';
    else
        engine = 'pa_wavplay';
    end
end
if ~exist('pageSize', 'var')
    pageSize = 0;
end
if pageSize > size(x,1)
    pageSize = 0;
end

% ------------------------------------------------ 2. determine in/out mode
if any(x(:)) && isempty(recordCH)
    % playback only
    ioMode = 1;
elseif ~any(x(:))
    % record only
    ioMode = 2;
else
    % playback and record
    ioMode = 3;
end

% ----------------------- 3. setup the devices if not specified by the user

% call AK setup functions
if isempty(playDevice) || (isempty(recordDevice) && ~isempty(recordCH))
    if strcmpi(engine, 'playrec')
        [playDevice, recordDevice, pageSize, fs, playCHmax, recCHmax] = AKplayrecSetup([], [], [], fs);
    elseif strcmpi(engine, 'pa_wavplay')
        [playDevice, recordDevice, playCHmax, recCHmax] = AKpa_wavplaySetup([], []);
    else
        error('AKio:Input', 'engine must be ''playrec'' or ''pa_wavplay''')
    end
end

% ------------------------- 4. default playback/record channels, and x size
% default channel
if (ioMode == 1 || ioMode == 3) && isempty(playCH)
    playCH = 1:size(x,2);
end
if (ioMode == 2 || ioMode == 3) && isempty(recordCH)
    recordCH = 1;
end

% set default maximum channels
if ~exist('playCHmax', 'var') 
    playCHmax = max(playCH);
elseif (ioMode == 1 || ioMode == 3) && isempty(playCHmax)
    playCHmax = max(playCH);
end
if ~exist('recCHmax', 'var')
    recCHmax = max(recordCH);
elseif (ioMode == 2 || ioMode == 3) && isempty(recCHmax)
    recCHmax = max(recordCH);
end

% check dimension of input buffer x
if size(x,2) < numel(playCH)
    x = repmat(x(:,1), [1 numel(playCH)]);
elseif size(x,2) > numel(playCH) && (ioMode ~= 2)
    x = x(:, 1:numel(playCH));
    warning('AKio:Setup', ['x has more channels than playCH. Only the first ' num2str(numel(playCH)) ' channel(s) are used'])
end

% ------------------------------------------------------ 5. check the input
% check if we have the needed devices
if ioMode == 1 && playDevice == -1
    error('AKio:Setup', 'The playDevice was not set - playback stopped')
elseif ioMode == 2 && recordDevice == -1
    error('AKio:Setup', 'The recordDevice was not set - recording stopped')
elseif ioMode == 3 && (playDevice == -1 || recordDevice == -1)
    error('AKio:Setup', 'The playDevice and or recordDevice were not set - playback and recording stopped')
end

% check if the desired channels are available
if (ioMode == 1 || ioMode == 3) && any(playCH > playCHmax)
    error('AKio:Setup', 'At least one playCH is higher than the maximum available channel number')
elseif (ioMode == 2 || ioMode == 3) && any(recordCH > recCHmax)
    error('AKio:Setup', 'At least one recordCH is higher than the maximum available channel number')
end
            
% -------------------------------------- 6. initialize playrec if necessary
if strcmpi(engine, 'playrec')
    %  reset playrec if neccessary
    if playrec('isInitialised')
        playrec('reset')
    end
    %  initialize
    if ioMode == 1
        playrec('init', fs, playDevice, -1, playCHmax, recCHmax, pageSize)
    elseif ioMode == 2
        playrec('init', fs, -1, recordDevice, playCHmax, recCHmax, pageSize)
    else
        playrec('init', fs, playDevice, recordDevice, playCHmax, recCHmax, pageSize)
    end
end


% ------------------------------------------- 7. audio playback / recording
if strcmpi(engine, 'playrec')
    
    % use playreac for audio IO
    if ioMode == 1
        AKplayrec('play', x, playCH);
    elseif ioMode == 2
        recording = AKplayrec('rec', x, playCH, recordCH, size(x,1));
    else
        recording = AKplayrec('playrec', x, playCH, recordCH);
    end
    
else
    
    % get the playback buffer
    paBuffer           = zeros(size(x,1), max(playCH));
    paBuffer(:,playCH) = x;
    
    % use pa_wavplay for audio IO
    if ioMode == 1
        pa_wavplay(paBuffer, fs, playDevice);
    elseif ioMode == 2
        recording = pa_wavrecord(min(recordCH), max(recordCH), size(paBuffer,1), fs, recordDevice);
    else
        recording = pa_wavplayrecord(paBuffer, playDevice, fs, 0, min(recordCH), max(recordCH), recordDevice);
    end
    
    % discard undesired recorded channels (pa_wavplay can only record
    % neighboring channels)
    if (ioMode == 2 || ioMode == 3) && size(recording,2) ~= numel(recordCH)
        recording = recording(:,recordCH-min(recordCH+1));
    end
end         

if ~exist('recording', 'var') && nargout
    recording = [];
end