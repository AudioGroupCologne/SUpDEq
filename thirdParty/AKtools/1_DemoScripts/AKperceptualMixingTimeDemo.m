% This matlab script accompanies the article:
%
% Lindau, A.; Kosanke, L.; Weinzierl, S.: 
% "Perceptual Evaluation of Model- and Signal-Based Predictors of the 
% Mixing Time in Binaural Room Impulse Responses", In:  J. A. Eng. Soc.,
% 60(11):887-898, (2012).
%
% It serves for computing estimates of the perceptual mixing time, i.e. 
% the instant a room's sound field is percieved as being diffuse. The 
% perceptual mixing time is estimated using model- and/or signal-based 
% predictors. Please note that - assuming a worst-case situation - 
% predictors were developed for the application in rectangular (i.e. 
% shoe-box-shaped) rooms. Their validity might thus be limited to such 
% rooms.
%
% As model-based predictors the enclosure's ratio of volume over surface 
% area V/S, or its volume V are used, respectively. If an impulse response 
% is available, the data-based predictor from Abel & Huang 
% (2006, criterion I, cf. above cited article) is be applied.
%
% In both cases, the perceptual mixing time is predicted as the threshold 
% value for
% a) the average assumed normal distributed population, and as a more
%    strict criterion,
% b) for the 95%-point of the assumed normal distributed population.
%
% In case of data-based estimation, the provided impulse responses should 
% fulfill three conditions. First they should be delivered in a 
% matlab-readable *.wav-format (single- or multichannel). Second, the 
% direct sound in the impulse response should be the strongest signal 
% component as the calculation conducts a peak detection with according 
% assupmtion. Third, per default the peak-detect algorithm assumes a 
% peak-to-noise-ratio of at least 30 dB (value can be changed in the script
% via the variable onset_threshold_dB).
%
% After running the code you will be guided through it in step-by-step 
% instructions.
%
% This code has been tested with Matlab vs. R2016a. Figure's format might
% be corrupted when submitting impulse responses with more than 2 channels.
% 
% Additional references:
% 
% Abel & Huang (2006): "A simple, robust measure of reverberation echo
% density", In: Proc. of the 121st AES Convention, San Francisco
%
% All impulse responses used in the study by Lindau and Kosanke can be
% downloaded from:
% http://www.ak.tu-berlin.de/menue/digitale_ressourcen/research_tools/mixing_time_prediction/
%
%
% A. Lindau, L. Kosanke, 2011
% alexander.lindau@tu-berlin.de
% audio communication group
% Technical University of Berlin

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
close all; clear;  clc

%% ----------------------------------------------- 1. processing parameters

N                  = 1024;  % window length (see Abel & Huang, 2006)
onset_threshold_dB = -30;   % peak criterion(onset_threshold*maximum)
peak_secure_margin = 100;   % number of samples used to protect peak from being affected, when elimiating time of flight


%% ----------------------------------- 2. model-based predictors: V/S and V
disp('----------------------------------')
disp('MODEL-BASED MIXING TIME PREDICTION')
disp('----------------------------------')
disp(' ')
disp('ROOM PROPERTIES (ENTER to skip model-based prediction)')

h = input('height in [m]: ');

if ~(isempty(h)) % ... for skipping...
    
    doModelBased = true;
    
    if  h<0 || ~(isnumeric(h))
        error('Insert positive, numeric values only.')
    end
    l = input('length in [m]: ');
    if  l<0 || ~(isnumeric(l))
        error('Insert positive, numeric values only.')
    end
    b = input('width in [m]: ');
    if  b<0 || ~(isnumeric(b))
        error('Insert positive, numeric values only.')
    end
    
    % compute perceptual mixing time (model_based)
    [tmp50_model_based, tmp95_model_based] = AKmodelBasedMixingTime(h,b,l);
    
else
    
    doModelBased = false;
    
    disp('Skipped model-based mixing time prediction.')
    tmp50_model_based = [];
    tmp95_model_based  = [];
end

%% -------------- 3. data based predictor: Abel & Huang (2006, criterion I)
disp(' ')
disp('---------------------------------')
disp('DATA-BASED MIXING TIME PREDICTION')
disp('---------------------------------')
disp(' ')

% select an impulse response
% - by default, we use an IR from the AKtools demo data
% - more IRs can be found at: http://www.ak.tu-berlin.de/menue/digitale_ressourcen/research_tools/mixing_time_prediction/
default_IR = fullfile(fileparts(which('AKtoolsStart.m')), '0_DemoData', 'BRIR FABIAN Teldex V3600 RT1_8s.wav');

[IR_name, IR_path] = uigetfile({'*.wav'; '*.ogg'; '*.flac'; '*.au'; '*.aiff'; '*.aifc'}, 'Select impulse response (press ''cancel'' to skip)', default_IR);

% check input
if ~IR_name
    
    doDataBased = false;
    
    display('Skipped data-based mixing time prediction.')
    tmp50_data_based        = [];
    tmp95_data_based        = [];
    tmp50_interchannel_mean_data_based   = [];
    tmp95_interchannel_mean_data_based   = [];
    echo_dens               = [];
    fs                      = [];
else
    
    doDataBased = true;
    
    % load IR
    IR_name = fullfile(IR_path, IR_name);
    [IR,fs] = audioread(IR_name);
    
    disp('IMPULSE RESPONSE PROPERTIES:')
    channel_number  = input('Channel number (ENTER to use all): ');
        
    % specific channel number chosen?
    if size(channel_number,1) == 1
        % check channel number
        if size(IR,2) < channel_number
            disp('Channel does not exist, channel 1 is used')
            channel_number = 1;
        end
        % use the chosen channel of the IR
        IR = IR(:,channel_number);
    end
    
    % ask for stopping time
    fprintf('\n')
    display('To increase calculation speed: Provide a maximum length in [ms]')
    display('of the impulse response to be taken into account (e.g. twice')
    display('the expected mixing time.)')
    display('ENTER to choose default (300 ms),')
    display('insert [0] for complete IR, or')
    stop_time = input('[-1] to stop after finding the mixing time (faster): ');
    if stop_time == 0
        fprintf('\nUsing complete impulse response for calculation (might take a while).\n \n')
        speedUp = false;
    elseif stop_time == -1
        fprintf('\nAbort calculation after finding the mixing time.\n \n')
        stop_time = 0;
        speedUp   = true;
    else
        if isempty(stop_time)
            fprintf('\n')
            display('Using 300 ms (default) as stopping time!')
            stop_time = 300;
            speedUp   = false;
        end
    end
        
    % compare IR-length and provided stop_time
    if stop_time ~= 0
         stop_time_samples = floor(stop_time/1000*fs);
         if length(IR) < stop_time_samples 
             fprintf('\nProvided impulse response is shorter than demanded stopping time.')
             display('Using complete impulse response for calculation (might take a while).')
             stop_time = 0;
         end
    end
        
    % cut IR (use IR only from onset position to stop_time) 
    IR = AKcutIRmixingTime(IR,fs,stop_time,onset_threshold_dB,peak_secure_margin);

    % compute perceptual mixing time (data_based)
    [tmp50_data_based, tmp95_data_based, tmp50_interchannel_mean_data_based, tmp95_interchannel_mean_data_based,echo_dens] = AKdataBasedMixingTime(IR,N,fs,peak_secure_margin,speedUp);

end

%% -------------------- 4. plot and save results if anything was calculated
if doModelBased || doDataBased
    
    disp(' ')
    disp('---------------------------------')
    disp('      Plot and save results      ')
    disp('---------------------------------')
    disp(' ')
    
    % save results as picture?
    do_print    = input('Save results as tiff to current directory [y/n] (ENTER = No)?','s');
    if strcmpi(do_print, 'y')
        do_print = 1;
    else
        do_print = 0;
    end 
    
    % run the plot function
    AKplotMixingTime(tmp50_model_based, tmp95_model_based, tmp50_data_based, tmp95_data_based, tmp50_interchannel_mean_data_based, tmp95_interchannel_mean_data_based, echo_dens,fs,do_print)
end