% Demonstration of (multi-channel) transfer function inversion using
% frequency dependent regularization. The compensated system - i.e. the
% convolution of the transfer function with the compensation filter - can
% be designed to approach either minimum or linear phase.
%
% This implementation of the regulated inversion has quite some parameters.
% It might be easiest to run the script and look at the reuslts to
% understand whats happening...
%
% In this example we will use regulated inverions to design inverse filters
% for headphone transfer functions. In this case the regularization is used
% to limit the gain of the compensation filters in order not to compensate
% sharp and deep notches (valleys) in the headphone transfer functions.
% However there are more areas for applying regularization. For instance
% when designing cross-talk filters for transaural binaural synthesis or
% when compensating the transfer function of loudspeakers in rooms.
%
% Sections 1-5 set the paramters for the inverions, section 6 runs it:
% 1. input/output data
% 2. input data pre-processing
% 3. target function
% 4. inversion / regularization parameters
% 5. plots
% 6. compute compensation filter
% 7. save intermediate results and documentation
%
% Zora Scharer , Alexander Lindau (2009):"Evaluation of Equalization
% Methods for Binaural Signals." 126th AES Convention, Munich, Germany.
%
% Alexander Lindau, Fabian Brinkmann (2012):"Perceptual evaluation of
% headphone compensation in binaural synthesis based on non-individual
% recordings." J. Audio Eng. Soc., 60(1/2), 54-62.
%
%
% F. Brinkmann, A. Lindau, Z. Schaerer
% TU Berlin, Audio Communication Group. 2008-2016

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

%% ---------------------------- 0a. select directory for saving the results

% (sub-directories will be made automatically)
% Note: you can replace uigetdir with a fixed path if you like
disp('----- Select a directory for saving the results -----')
s.save_dir  = uigetdir(fullfile(fileparts(which('AKtoolsStart'))), 'Select a directory for saving the results');

%% --------------------------------------- 0b. give the path to this script
%  this is done because the script will be copied for documentation

% Note: you can replace uigetfile with a fixed path if you like
disp('----- Point to this script for saving it for documentation -----')
[script_name, script_dir]  = uigetfile(fullfile(fileparts(which('AKtoolsStart'))), 'Point the the script you are using');
s.script_dir = fullfile(script_dir, script_name);

clear script_dir script_name

%% --------------------------------------------------- 1. input/output data
%  Define your input data and location and names for your output data in
%  this section

% load input data
% - we use headphone impulse responses (HpIRs) for this purpose
% - you can also delete the next three lines and read wav files using in s.data_dir
h    = SOFAload('HpIRs.sofa');
s.fs = h.Data.SamplingRate;     % set the sampling rate
h    = shiftdim(h.Data.IR, 2);  % format the HpIRs (12 HpIRs for each ear, with a length of 4096 samples)
    

% directory containing the raw transfer functions to be inverted.
% All wave-files in this directory will be read!
% This is ignored if you manually pass the IRs, as we do it in this demo
s.data_dir   = '..';

% save results 
s.save_filter = 2;              % save compensation-filter, 0: don't save, 1: fWonder-Format, 2: multichannel wav-file, 3: multichannel mat-file, 4: SOFA 'GeneralFIR' convention, 5: SOFA 'SimpleFrieFieldHRIR' convention for use with a modified version of the SSR
s.save_name   = 'HpIR_filter';  % name for saving compensation-filter (do not append 'wav' or 'mat')
s.nbits       = 32;             % word length used for wavfile saving (Note: use 32 bit if you save to wav files, to avoid clipping)
s.att         = 0;              % optional attenuation of final compensation filter by -s.att dB before saving
s.fade_in     = 64;             % duration of sin^2 fade in on compensation filter in samples
s.fade_out    = 64;             % duration of cos^2 fade out on compensation filter in samples

% documentation
s.save_doc    = 1;              % copy filter calculation scripts, and intermediate results to s.save_dir


%% ------------------------------------------- 2. input data pre-processing
%  Specifiy the pre-processing steps of the input data in this section

s.Nsamples = 2^12;              % input data will be truncated to s.Nsamples. Compensation filters will have the same length.
s.win_len  = 256;               % fade out duration in samples after shortening input data to s.Nsamples
s.dyn      = false;             % limit the dynamic range of the input data and target function to -abs(s.dyn) dB FS.
                        		% might be helpful in case the filter design is disturbed, but won't be needed in most cases
                                
s.avg_type     = '_complex';    % If multiple measurements are passed per channel, the input data will be averaged
                                % '_complex' - average complex spectra
                                % '_time'    - average impulse responses
s.time_align   = false;         % time align impulse responses before averaging based on time-of-arrival (TOA) estimation (e.g. to compensate a varying sound-card latency)
s.threshold    = -20;           % threshold in dB that used for detecting the TOA in the input data.
                                % The first sample that exceeds s.threshold defines the TOA
s.peak_protect = 20;            % if s.time_align is set to true, the onsets in the input data will be shifted to s.peak_protekt samples.
                                % You should keep some samples to avoid cutting the beginning of the impulse responses

s.APPLY_SMOOTHING = false;      % fractional octave smoothing of  averaged transfer function before inversion. Phase is made minimum phase in this case.
s.frac            = 6;          % fraction for octave smoothing (e.g. 3 -> third octave smoothing).
                                % Attention: Value will effect final filter's latency ! Check smallest possible values carefully, use do_plot = 1


% The level of the averaged input data will be normalized.
% This is important because the regularization will only be active, if the
% magnitude response of the transfer function is below the target function
% (see below)
% Interchannel differences (e.g due to different sensitivities of left and
% right headphone) are saved and applied to compensation filters
% afterwards. Normalization should be done at a frequency where the
% transfer functions show a constant magnitude response.
s.f_norm       = 200;              % the level at the frequency s.f_norm in Hz will be used for normalization
s.f_norm_frac  = 1;                % width in fractional octaves around s.f_norm that will be used for normalization
                                   % (.5= 2 octaves, 3=third octave, etc. Range is shown in plot2)
s.g_norm       = 0;                % level in dB to which the average transfer function is normalized


% filter to be applyied before and after inversion. the pre-filter
% compensates the transfer function of the measurement situation
% (microphones, pre-amps, etc.), the post-filter that of the playback
% situation (microphones included in binaural recordings, etc.).
% The filters have to be saved in a mat-file that contains a variable
% called pre_post_filter. An example for that can be found in the
% 0_DemoData folder ('MicrophoneFilter.mat')
% - false if nothing should be applied
% - full path, or filename,
s.prefilter  = false;
s.postfilter = false;           % if the postfilter is applied, the length of the compensation filter will be s.Nsamples + length(postfilter) - 1


%% ----------------------------------------------------- 3. target function
%  The convolution of the input data with the compensation filter will
%  apprach the target function that is defined in this section. In this
%  demo we will use either a band-pass or a constant spectrum but different
%  targets (e.g. a diffuse field response) can be implemented quite easily.
%  
%  In many cases you might want to restric the band-width of the inversion
%  process to a range where the input data is valid. This can either be
%  done by defining a band pass target function or by applying
%  regularization at high and low frequencies.

% define target function to which the compensation filter is fitted
% (usefull to restrict the range within the filter should work)
% - define butterworth band-pass by 4-element vector
%   (for high-pass/low-pass only, pass 0 for hp_order/lp_order)
%   [hp-order, hp_cut-off_freq, lp-order, lp_cut-off_freq]
% - do not apply bandpass by passing
%   false
% - specifiy BKamp preset by string identifier:
%   (this is only interesting if you use the extra aural BK headphones from
%   the audio communication group)
%   'fr' (fullrange)
%   'sub50'
%   'sub85'
%   'sub120'
%   'sub150'
s.target      = [2 30 1 16000]; % for this demo use the target for all methods just for simplicity
s.phase_type  = 'min_phase';    % phase of target function: min_phase, lin_phase
s.Nfft_double = 6;              % zero-pad before generating minimum-phase to make it more robust. 1=double length, 2=quadruple lenght, etc.
s.margin      = 1200;           % The compensation results will approach the phase behaviour defined by s.phase_type.
                                % If you chose 'min_phase' the impulse response of the target function is circshifted to s.margin samples.
                                % The compensation filter might be compromised if this value is to low.
                                % 1200 samples will work in almost all cases, smaller values often work as well.
                                % Try to reduce this value to for application that demand low latency!


%% ------------------------------- 4. inversion / regularization parameters
% choose method for inversion /regularization and set parameter in the
% sections below.
%
% If you choose a method that uses regularization it will roughly work like
% this:
% You define a regularization function. The gain of the compensation filter
% is limited at frequencies where the regularization function has large
% vales, i.e. 0 dB FS. The regularization weight beta will determine how
% stringly the inversion is limited. In most cases, the regularization
% weights should be somwhere between 0.5 ... 0.003
%
% The following methods are implemented in AKtools:
%
% 1 - regularized least mean squares (LMS) inversion with a high-shelf for
%     regularization
% 2 - regularized LMS inverions using the inverted, and smoothed input data
%     for regularization
% 3 - LMS inversion of the smoothed input data (no regularization)
% 4 - regularized LMS inversion using finely and coarsly smoothed versions
%     of the input data for regularization. This regularization function
%     will have high values at the frequencies where the input data shows
%     sharp and deep notches (valleys). This is similar to compare and
%     squeese from Swen Mueller
% 5 - direct inversion of the fractional octave smoothed input data with
%     limited inversion dynamic
% 6 - regularized LMS inversion with a regularization function generated
%     from parametric equalizers, highshelfs, and low-shelfs
%
% At least for headphones methods 4 and 6 tend to give the best results.

s.inversion_method        = 6;

switch s.inversion_method
    case 1
        s.reg_shelve_gain  = [15 15];      % channel wise shelving filter gain
        s.reg_shelve_freq  = [4000 4000];  % channel wise shelving filter freq. @ half gain
        s.reg_beta         = [.3 .3];      % channel wise regularization weights
    case 2
        s.reg_frac         = 6;            % fraction for octave smoothing (e.g. 3 -> third octave smoothing)
        s.reg_beta         = [.5 .5];      % channel wise regularization weights
        s.reg_align        = [300 16000];  % set regularization function flat below s.inv_data.align(1) and above align(2) Hz
        s.reg_inv_dyn      = 30;           % constrain inversion dynamic of regularization function in dB

    case 3
        s.frac             = 3;            % fraction for octave smoothing (e.g. 3 -> third octave smoothing)
                                           % s.APPLY_SMOOTHING and s.frac will be overwritten)

    case 4
        s.reg_frac_fine    = 3;            % fraction for fine octave smoothing
        s.reg_frac_rough   = 1;            % fraction for coarse octave smoothing
        s.reg_beta         = [.2 .2];      % channel wise regularization weights
        s.reg_align        = [300 20000];  % set frequency response flat below s.inv_data.align(1) and above align(2) Hz

    case 5
        s.inv_data.frac    = 24;           % fraction for octave smoothing (e.g. 3 -> third octave smoothing)
                                           % s.APPLY_SMOOTHING and s.frac will be overwritten)
        s.inv_data.inv_dyn = 30;           % constrain inversion dynamic in dB

    case 6
        s.reg_beta         = [.3 .3];      % channel wise regularization weights
        s.reg_dyn          = [10; 10];     % channel wise dynamic of the regularization function
                                           % (dynamic is limited by clipping the regularization function)
                                           % By clipping a parametric equlizer the regularizatoin function
                                           % can have plateaus of equal regularization which might be
                                           % helpful in case of broader notches
        
        % specify arbitrary number of parametric equalizers (PEQ) and self
        % filters by:
        % 1. Channel to be applied to (referring to the input data)
        % 2. Type: 'PEQ', 'LS' (low-shelf), or 'HS' (high-shelf)
        % 3. PEQ center frequency, or shelf characteristic frequency (see 5.)
        % 4. Gain
        % 5. PEQ Quality
        %    Shelf characterisitc frequency type
        %    - NaN   :3 specifies the shelf's mid gain frequency
        %    - '3dB' :3 specifies the shelf's 3dB cut-off frequency
        s.reg_data    = {
                          % 1 'LS'  40    10 NaN;...  % <- we don't need the shelf filter if we apply a band pass target function. But you can go ahead and change that...
                          1 'PEQ' 8323  20 10 ; ...
                          % 1 'HS'  13000 7  '3db';...
                          % 2 'LS'  40    10 NaN;...
                          2 'PEQ' 8710  20 10 ; ...
                          2 'PEQ' 12870 20 10 ; ...
                          % 2 'HS'  13000 7  '3db';...
                         };        
end

%% --------------------------------------------------------------- 5. plots
%  different plots can be selected to find the otpimal parameters for the
%  inversion
%
%  You should consider the following things (they are at least important
%  for headphone compensation)
%  - Your pre-processed input data shown in plots 1, and 2 should look
%    clean (e.g. all impulse responses should start at the same time, there
%    should be no extreme outliers)
%  - Your regularization function shown in plot 3 (if used) should fit your
%    data, i.e. have large values only at frequencies where you want to
%    restrict the inversion
%  - your target function shown in plot 4 (if used) should correspond to
%    the bandwidth of your input data.
%  - The compensation results applied to the averaged input data as shown
%    in plot 5 should approach the target function
%  - The inverse filter shown in plot 6 should be free of peaks with high
%    gains. Be also careful with boosting low frequencies.
%  - The resuts shown in plot 7 should be as close to 0 dB FS as possible.
%
%  Trouble shooting
%  - if your compensation filter has a strong low frequency ripple either
%    increase s.margin, apply regularization at low frequencies (e.g. by
%    using s.inversion_method=6), or apply a band-pass target function
%  - if your compensation results looks low passed, try smaller values for
%    regularization weight beta, or apply a strong frequency dependant
%    regularization function (e.g. s.inversion_method=4 or 6)
%  - if your compensation filter exhibits high and small band peaks
%    increase the regularization weights, or change the regularization
%    function

% plotting and printing
s.plot_shorten        = 1;      % plot 1: results of toa alignement, and shortening the input data
s.plot_avg            = 1;      % plot 2: results of averaging the input data
s.plot_regularization = 1;      % plot 3: regularization function
s.plot_target         = 1;      % plot 4: target function
s.plot_compensation1  = 1;      % plot 5: summary of inversion process
s.plot_inverseFilter  = 1;      % plot 6: inverse filter
s.plot_compensation2  = 1;      % plot 7: compensation filter applied to input data
s.do_print            = 1;      % save plots as pdf's to s.save_dir

%% ----------------------------------------- 6. compute compensation filter
% calculate headphone compensation filter
if exist('h', 'var')
    data = AKregulatedInversion(s, h);
else
    data = AKregulatedInversion(s);
end


%% ------------------------- 7. save intermediate results and documentation
% copy matlab files used for filter generation and save some preliminary
% results
if s.save_doc
    % make filter design documentation directory
    if exist(fullfile(s.save_dir,'filterdesign_doc'), 'dir') ~= 7
        mkdir(fullfile(s.save_dir,'filterdesign_doc'))
    end
    % copy scripts
    copyfile(s.script_dir, fullfile(s.save_dir,'filterdesign_doc','BKP_AKregulatedInversionDemo.m'))
    copyfile(which('AKregulatedInversion'), fullfile(s.save_dir,'filterdesign_doc','BKP_AKregulatedInversion.m'))
    copyfile(which('AKffdInverseFilter'), fullfile(s.save_dir,'filterdesign_doc','BKP_AKffdInverseFilter.m'))
    % remove large data from the save data
    if isfield(data, 'raw')
        data = rmfield(data, 'raw');
    end
    if isfield(data, 'prefilter')
        data = rmfield(data, 'prefilter');
    end
    if isfield(data, 'aligned')
        data = rmfield(data, 'aligned');
    end
    % save data
    save(fullfile(s.save_dir,'filterdesign_doc','data'), 'data')
end
