% Generic script for simple impulse response (IR) measurement using a swept
% sine excitation signal.
%
% You need a sound card with an ASIO4all (http://www.asio4all.com) driver
% if you run this in windows!
%
% - For the first measurement you will have to evaluate this script block
%   by block (cmd+enter on mac; ctrl+enter on linux/windows)
% - After that you can evaluate single blocks to only change certain
%   parameters. For instance you can execute Block 7 multiple times to
%   measure multiple impulse responses
%
%  1. data/global settings
%  2. initialize audio interface
%  3. generate excitation signal
%  4. setup post processing
%  5. reference measurement
%  6. level calibration
%  7. impulse response measurement
%
% v1.0 2016/05 inital development
%              fabian.brinkmann@tu-berlin.de, david.ackermann@tu-berlin.de
% v1.1 2017/11 increased usablilty
%              helmholz@campus.tu-berlin.de, fabian.brinkmann@tu-berlin.de

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
close all; clear; clc

%% ------------------------------------------------ 1. data/global settings
% sampling rate in Hz
s.fs = 44100;

% select directory for saving data
data.dir    = [];       % full path of the directory for saving the data (if this is empty the directory is selected in a GUI)

% select which data to save
data.IR     = true;     % save IR after deconvolution (variable 'ir')
data.raw    = false;    % save recorded sweep and excitation signal (variables 'raw' and 'sweep')
data.ref    = true;     % save reference measurement (variable 'reference')
data.calib  = true;     % save calibration measurement (variable 'calibration')
data.meta   = true;     % save meta data (variables 's', 'm', 'x', and 'data')
data.plot   = true;     % save plots
data.close  = true;     % close plots automatically

data = AKmeasureCheckDir(data);


%% ------------------------------------------ 2. initialize audio interface
% select audio engine based on operating system
if isunix
    s.audio_engine    = 'playrec';
else
    s.audio_engine    = 'pa_wavplay';
end

AKmeasureIO


%% ------------------------------------------ 3. generate excitation signal
% settings for excitiation and deconvolution
% If you want to applay different/additional settings you will also have to
% change the call for AKdevonc in this script and in AKmeasureIR.m
x.NFFT     = 18;             % order (duration in samples = 2^x.NFFT)
x.highPass = [20 8 1];       % lower cut-off frequency given by pass-band and stop-band of a low-pass [Hz]
x.lowPass  = false;          % higher cut-off frequency is deactivated by default
x.start    = 60;             % silence before start of sweep in ms
x.end      = 2000;           % silence after end of sweep in ms (must be long enough to allow for the decay of the room)
x.dynamic  = 40;             % dynamic range [dB] of the convolution. Dynamic of inverted excitatin signal is restricted accordingly


% generate and plot excitation signal
s.sweep = AKsweepFD('NFFT', x.NFFT, 'fs', s.fs, 'preset', 'log', ...     % duration, and type
                    't_start', x.start, 't_gap', x.end, ...              % silence before and after sweep in ms
                    'highPass', x.highPass, 'lowPass', x.lowPass, ...    % restrict band width
                     'do_plot', true);

if data.plot
    print('-dpdf', fullfile(data.dir, 'Plots', 'excitiationSignal'))
end

% invert and plot excitation signal
[~] = AKdeconv(s.sweep, s.sweep, s.fs, 'x_inv_dyn', x.dynamic, 'plot_inv', 1);
if data.plot
	print('-dpdf', fullfile(data.dir, 'Plots', 'excitiationSignalInverted'))
end

% close plots
if data.close
    close all
end

%% ----------------------------------------------- 4. setup post processing
%  Sometomes it is nice to have some automated post-processing. This can be
%  found in AKmeasureIR. You can implement more stuff here, if you want...

s.subSonic = true;  % apply 20 Hz high pass to filter out noise
s.N        = false; % shorten to p.N samples before saving IRs, false for keeping original length

%% ----------------------------------------------- 5. reference measurement
% cancel out some effects of your measurement chain by selecting one of the
% following referenceTypes. The reference will always contain a single
% channel.
%  false     - nothing is done (if desired set s.referenceLatency below manually)
%  'latency' - the latency is measured and corrected using a circular shift
%  'complex' - the complex frequency response is measured and corrected
s.referenceType        = false;  

% output to input latency in samples.
% This can be set if you want to correct the latency with out measuring
% (set s.referecyType = false in this case. s.referenceLatency will be
% overwritten if s.referenceType = 'latency')
s.referenceLatency     = 0; 

s.referenceChOut       = 2;     % output channel for reference measurement
s.referenceChIn        = 2;     % input channel for reference measuremnt
s.referenceLevelOut_dB = -10;   % output level in dB FS for reference measurement    

if s.referenceType
    AKmeasureReference % see inside the script for more information
end


%% --------------------------------------------------- 6. level calibration
% Calibrates the level of the input chain. If calibration is applied a
% value of 1 in the impulse response will refer to 1 Pascal.

% Set the calibration mode:
%  false     - do not calibrate the input chain
%  true      - apply calibration with a test tone (e.g. a calibrator)
%  'numeric' - manual calibration by passing numeric value that specifies 
%              the dB FS (peak not RMS!) value that is expected at an input signal with an
%              amplitude of 1 Pascal (can be calculated from the
%              microphone, pre-amp, and sound-card specifications)
s.calibrate             = false;    % true to apply calibration, false to omit
s.calibrateChIn         = 1;        % input channel for calibration
s.calibrationFrequency  = 1000;     % calibration frequency in Hz
s.calibrationLevel      = 94;       % calibration level in dB SPL

if s.calibrate
    AKmeasureCalibrate
end


%% ------------------------------------------------------ 7. IR measurement
%  This will start the measurement. Make sure the output level is ok, and
%  does not blow your ears!

% Define input and output channels
s.chOut = 1;       % list of output channels e.g. [1 3].
                   % Output channels can measured be measured one after
                   % another or at the same time -> see s.chMode below.
s.chIn  = 1;       % list of input channels e.g. [1 3].
                   % Input channels are recorded simultaneously!

% Define how the output channels are measured
s.chMode= 'single'; % 'single' each channel is measured separately
                    % 'all'    all channels are measured at the same time

% Level and averaging
s.levelOut_dB   = -40; % output level in dB FS
s.clipReduction = 3;   % automatic reduction of level in dB if clipping occured
                       % (if gain reduction occurs, it is compensated automatically after the measurement)
s.average       = 1;   % number of averages (doubling the averages increases the
                       % signal-to-noise ratio by 3 dB)

% deconvolution settings
s.deconvType = 'lin';  % 'lin': non-harmonic distortion products at the end of the IR are discarded
                       % 'cyc': non-harmonic distortion products are kept
                       
%  some information about the measurement itself. You can add anything you
%  want here. It will be saved seperately for each meesurement.
m.contact     = 'Christiane Schmidt';
m.room        = 'reverberation chamber TU Berlin';
m.temperature = 17;      % in degree celcius
m.humidity    = 47;      % relative humidity
m.soundcard   = 'RME Fireface UCX @ 4 dBU outputSensitivity, 4 dBU inputSensitivity, SN:112';
m.amplifier   = 'none';
m.micPreamp   = 'LP C360, 24 dB gain, SN:112';
m.source      = 'QSC K8, SN:112, gain @ max, EQ flat';
m.receiver    = 'NTI MA220, SN:112';

% run the measurement
AKmeasureIR        % see inside the script for more information

%  Save the data (as specified in section 1 data/global settings)
%  You will be asked to name the files each time. If files with the same
%  name already exist, you will be asked if you want to override them.
AKmeasureSaveData
