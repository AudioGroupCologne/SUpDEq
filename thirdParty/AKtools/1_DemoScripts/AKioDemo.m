% This script demonstrates audio playback and recording with AKtools using
% AKio.m which is short for AKin/out.
%
% You need a sound card with an ASIO4all (http://www.asio4all.com) driver
% if you run this in windows!
%
% For audio playback and recording external audio engines (pa_wavplay, and
% playrec) are used. Playrec can be used from Windows and Mac, hoever might
% produce glitches in Windows. pa_wavplay can be used from Windows, only.
%
% 0. warning - watch your audio playback level
% 1. audio playback
% 2. audio recording
% 3. audio playback and recording
%
% See AKmeasureDemo.m for a demo on how to measure impulse responses
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
close all; clear; clc
%#ok<*NASGU>

%% ------------------------------------------------------------- 0. warning
warning('This script will play back audio - watch your output level!!!')

%% ------------------------------------------------------ 1. audio playback

% In the simplest case we only pass an audio signal to AKio.m, which than
% selects the audio engine by itself and will show a GUI to set up audio
% devices. Audio is played back via the first channel in this case.

% generate a 1 kHz sine of 1 second duration with -40 db FS
x = AKsine(1, 1000) * 0.01;

% play it back
AKio(x)

%% We can also set different output channels, and play back multi-channel 
%  signals

% generate a 1 kHz sine of 1 second duration with -40 db FS
x = AKsine(1, [1000 1050]) * 0.01;

% play it back
AKio(x, [1 2])

%% To avoid to setting-up the audio devices via the GUI, we have to pass a
%  playDevice and recordDevice (see help AKio.m). We can also pass the
%  sampling rate and audio engine (playrec or pa_wavplay) if we want.
%  NOTE: this might not work if you dont have a play device with ID=1. You
%        can get the available devices from
%        AKplayrecSetup (mac default), and
%        AKpa_wavplaySetup (windows default)

AKio(x(:,1), 1, [], 44100, 1, [])

%% ----------------------------------------------------- 2. audio recording

% In the simplest case we only pass the number of samples that we want to
% record. Now the first channel is used.
recording = AKio(44100);

%% Again, we can specify the recordng channel(s)
recording = AKio(44100, [], [1 2]);

%% ---------------------------------------- 3. audio playback and recording

% For recording audio, we have to at least specify an audio signal and the
% playback and record channel(s)
recording = AKio(x(:,1), 1, 1);