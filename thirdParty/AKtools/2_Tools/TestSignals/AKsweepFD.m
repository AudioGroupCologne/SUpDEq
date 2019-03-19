% [sweep, group_delay, f_group_delay] = AKsweepFD(varargin)
%
% Frequenzy domain sine sweep generation [1]. The phase is calculated from
% user specified magnitude response. Time signal is obtained by inverse
% fourier transform. Spectral coloration can be achieved by specifying
% filters (low- and high-shelf) or passing noise and speaker impulse
% responses (for sweep measurements with constant SNR). Frequency range can
% be limited by a band-pass. Arguments are passed as 'Atribute', value
% pairs.
%
% EXAMPLES FOR CALLING:
% AKsweepFD                      all default parameters are used
% AKsweepFD('NFFT', 19)          FFT order = 19
% AKsweepFD('preset', 'log')     logarithmic sweep
% AKsweepFD('preset', 'perfect') perfect sweep
%
% See AKtestSignalsDemo.m for examples
% 
%
% I N P U T: (DEFAULT VALUES GIVEN IN BRACKETS):
% preset (none)         - For fast generation of typical sweep types. Note
%                         that some user/default parameters will be changed
%                         lin: linear sweep
%                         log: logarithmic sweep
%                         perfect: perfect sweep
%
% fs (44100)            - sampling frequency in Hz
% NFFT (18)             - length of output vector:
%                         if NFFT <= 25 the output vector has
%                         length 2^NFFT
%                         if NFFT >= 26 the output vector has lentgth NFFT
%                         and NFFT must be even!
% NFFT_double (1)       - double fft size during sweep calculation. More
%                         robust and thus recommended (0 or 1)
% t_start (250)         - time in ms before the sweep starts
% t_gap (1000)          - silence in ms that is inserted after the sweep to
%                         allow for high frequencies to decay
%
% highPass [20 8 1]     - restricts sweep band width with Butterworth
%                         high pass. Three element vector:
%                         [cut-off frequency in Hz, Order, M]
%                         The high pass is applied M times. Pass false, to
%                         disable.
% lowPass (false)       - restricts the sweep band width with a Butterworth
%                         low pass (c.f. highPass)
%
% lowshelf_apply (1)    - bass emphasis with 2nd order low-shelv, 
%                         Q = 1/sqrt(2). Lowshelv specs:
% lowshelf_f (500)      - mid-gain frequency in Hz;
% lowshelf_G (20)       - gain in dB
%
% highshelf_apply (0)   - high frequency emphasis with 2nd order high-shelv
%                         Q = 1/sqrt(2). Highhelv specs:
% highshelf_f (5000)    - mid-gain frequency in Hz;
% highshelf_G (-10)     - gain in dB
%
% noise_ir (0)          - single channel impulse response or absolute path
%                         to wav file of measured noise floor. Will be
%                         inverted, and applied to the sweep magnitude
%                         response. Length is adjusted to match NFFT and
%                         NFFT_double
% noise_dynamic (20)    - values that are more than noise_dynamic dB below
%                         maximum of noise freq. resp. will be clipped.
%                         frequency response is evaluated in limits given
%                         by highPass, and lowPass
% noise_smoothing(3)    - fraction of octave smoothing that is applied to
%                         noise magnitude response
%
% speaker_ir (0)        - single channel impulse response or absolute path
%                         to wav file of measured loudspeaker data. Will be
%                         applied by means of emphasis (after inversion).
%                         Length is adjusted to match NFFT ans NFFT_double
% speaker_dynamic (20)  - values that are more than speaker_dynamic dB
%                         below maximum of speaker freq. resp. will be
%                         clipped (before inversion). frequency response is
%                         evaluated in limits given by highPass, and
%                         lowPass
% speaker_smoothing(3)  - fraction of octave smoothing that is applied to
%                         speaker magnitude response
%
% env_manipulation (0)  - sometimes the sweep does not have a constant
%                         envelope, but overshoots at start and ending.
%                         This forces the envelope to be constant but
%                         changes the magnitude spectrum. Obacht - Harter
%                         Fusch!
%
% wav_save (0)          - string specifying filename with absolute path
%                         (with or without '.wav', 0: not saving as wav)
% wav_Nbit (32)         - bit depth of wav file
%
% log_save (0)          - string specifying filename with absolute path
%                         (with or without '.txt', 0: not saving, IF SWEEP
%                          IS SAVED AS WAV, LOG IS WRITTEN TO SAME
%                          DIRECTORY WITH SAME NAME BY DEFUALT)
%
% do_plot (1)           - plot time signal, magnitude response, and
%                         time until the end of the sweep in dependency of
%                         the frequency
%                         (MUST HAVE while designing a sweep!!!)
%                         1             : print to screen 
%                         absolute path : is saved as pdf file
%                         (IF SWEEP IS SAVED AS WAV, PLOTs ARE BY DEFUALT
%                         SAVED TO SAME DIRECTORY AND IDENTICAL NAME)    
% do_advanced_plot (0)  - Plots group delay of non-harmonic distortion
%                         products (0, 1)
% verbose (0)           - prints the processing steps to the screen (0, 1)
%
%
% O U T P U T:
% sweep                 - time signal
% group_delay           - surprise - its the group delay (one-sided)
% f_group_delay         - frequency vector for group delay
%
%
% [1] Swen Mueller, Paulo Massarani: "Transfer function measurement with
%     Sweeps. Directors's cut including previously unreleased material and
%     some corrections." J. Audio Eng. Soc. \emph{(Original release)}.
%     49(6):443-471, June, 2001.
%
% NOTE: Default parameters have only been testet for fs = [44100 48000]
%
% v.1 2011/10 F. Brinkmann, A. Giese, F. Schultz, A. Lindau,
%             Audio Communication Group, TU Berlin

% PROCESSING STEPS:
% 1.  parse input data
% 2.  setup_sweep type
% 3.  calculate helpful variables
% 4.  calculate magnitude response
% 5.  caluclate phase response
% 6.  calculate sweep in time domain
% 7.  envelope manipulation
% 8.  assign output parameters
% 9.  save wav and write log
% 10. plots
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

%-------------------------------------------------------------------------%
function [varargout] = AKsweepFD(varargin)

inputnames = ...
{'fs', 'NFFT','NFFT_double', 't_start', 't_gap',... 
 'highPass', 'lowPass', ...
 'lowshelf_apply', 'lowshelf_f', 'lowshelf_G',...
 'highshelf_apply', 'highshelf_f', 'highshelf_G',...
 'noise_ir', 'noise_dynamic', 'noise_smoothing',...
 'speaker_ir', 'speaker_dynamic', 'speaker_smoothing',...
 'env_manipulation',...
 'wav_save', 'wav_Nbit', 'log_save',...
 'do_plot', 'do_advanced_plot' , 'preset', 'verbose'};

% Define default values
def.fs          = 44100;        % sampling frequency in Hz
def.NFFT        = 18;           % fft order
def.NFFT_double = 1;            % double fft size during sweep calculation to reduce fold backs
def.t_start     = 250;
def.t_gap       = 1000;

def.highPass        = [20 8 2];
def.lowPass         = false;

def.lowshelf_apply  = 1;     % bass emphasis with 2nd order low-shelv
def.lowshelf_f      = 500;   % mid-gain frequency in Hz;
def.lowshelf_G      = 20;    % gain in dB

def.highshelf_apply = 0;     % high frequency emphasis with 2nd order high-shelv
def.highshelf_f     = 5000;  % mid-gain frequency in Hz;
def.highshelf_G     = -10;   % gain

def.noise_ir          = 0;    % apply moise floor by means of emphasis
def.noise_dynamic     = 20;   % maximum dynamic of noise floor
def.noise_smoothing   = 6;    % frequency smoothing

def.speaker_ir        = 0;    % apply moise floor by means of emphasis
def.speaker_dynamic   = 20;   % maximum dynamic of speaker floor
def.speaker_smoothing = 6;    % frequency smoothing

def.env_manipulation = 0;    % harter Fusch

def.wav_save        = 0;     % 0 or string specifying absolute path (with or without '.wav')
def.wav_Nbit        = 32;    % Bits of wav file

def.log_save        = 0;     % 0 or string specifying absolute path (with or without '.txt')

def.do_plot          = 1;      % plot time signal and magnitude response (MUST HAVE while designing a sweep!!!)
def.do_advanced_plot = 0;      % phase response (wrapped and unwrapped) and group delay, for better understanding
def.preset           = 'none'; % sweep type
def.verbose          = 0;


%% ---------------------------------------------------------- 1. parse input
% Check the number of arguments to be even
if mod(size(varargin,2),2)
    error('No even number of arguments. Check submitted attribute/value pairs.')
end

% Stock the values
for n=1:2:size(varargin,2)-1
    % checks for all possible attribute/value pairs
    is_parameter = 0;
    for m=1:size(inputnames,2) 
        % detect submitted attribute/value pairs
        if strcmp(inputnames{m},varargin{n})
            % create input-variables from submitted attribute-value-pairs (supress output)
            [~] = evalc([inputnames{m},'=','varargin{n+1}']);
            is_parameter = 1;
        end
    end
    if ~is_parameter
        error(['No such parameter: ' varargin{n}])
    end
end

% Create missing input-variables with default-values
for m=1:size(inputnames,2)
    default = ['def.' inputnames{m}];
    % if input-variable m hasn't already been defined...
    if ~exist(inputnames{m},'var')
        [~] = evalc([inputnames{m},'=',default]); % (supress output)
    end
end

% delete unnecassary variables to have a clean workspace inside the
% function (better for debugging)
clearvars def default in_args inputnames m n ninput trash varargin

%% --------------------------------------------------- 2. setup: sweep type
%TYPE OF SWEEP
if verbose
    fprintf('\n\n')
end

switch preset
    case 'none'
    case 'perfect'
        NFFT_double = 0;
        t_start = 0;
        t_gap = 0;
        highPass = false;
        lowPass  = false;
        lowshelf_apply = 0;
        highshelf_apply = 0;
        warning('AKseepFD:presets', 'Preset ''perfect'': Following Parameters were set to 0 (default/user input irgnored): NFFT_double, t_start, t_gap, highPass, lowPass, lowshelf_apply, highshelf_apply')
    case 'lin'
        lowshelf_apply = 0;
        highshelf_apply = 0;
        warning('AKseepFD:presets', 'Preset ''lin'': Following Parameters were set to 0 (default/user input irgnored): lowshelf_apply, highshelf_apply')
    case 'log'
        lowshelf_apply = 0;
        highshelf_apply = 0;
        warning('AKseepFD:presets', 'Preset ''log'': Following Parameters were set to 0 (default/user input irgnored): lowshelf_apply, highshelf_apply')
    otherwise
        error(['''' preset ''': No such preset'])
end
    

%% ----------------------------------------- 3. calculate helpful variables
if verbose
    fprintf('\n\n');
    disp('CREATING SWEEP WITH ARBITRARY AMPLITUDE SPECTRUM');
end
% Mueller + Massarani 2001, S.38:
%...a safe way to avoid contaminating the late decay of the tail by low-fre
%   quency components is to simply choose an FFT block length that is at
%   least double the desired sweeep length.
%
%   double NFFT-Size: abrupt sweep will cause less folding backs 
%   into the tail of the sweep
if NFFT_double
    if NFFT <= 25 %#ok<*UNRCH>
        n=NFFT+1;
    else
        n = 2*NFFT;
    end
else
    if NFFT <= 25
        n=NFFT;
    else
        if mod(NFFT,2)
            error('If NFFT is larger than 25, it must be even')
        else
            n = NFFT;
        end
    end
end

if (NFFT <= 26 && NFFT_double) || NFFT <= 25
    N = 2^n;      % number of samples
else
    N = n;
end
df  = fs/N;     % spacing between frequency bins of FFT
nyq = N/2+1;    % index of nyquist frequency (Only holds true vor even N!!!)
T   = N/fs;     % end time in seconds

%% ---------------------------------------- 4. calculate magnitude response
if verbose
    disp('calculating magnitude response...')
end

% APPLY COLORATION FOR LOG SWEEP
if strcmpi(preset, 'log')
    h_sweep = zeros(N/2,1);
    f_min = df;
    f_max = fs/2-df;
    h_sweep(2:N/2) = 1./sqrt(2*pi*(f_min:df:f_max));   
else
    h_sweep = ones(N/2,1);
end
clearvars f_min f_max 

% LOWSHELVING
if lowshelf_apply
    if verbose
        disp('low shelving...');
    end
    [b, a] = AKhighshelve2(lowshelf_f, fs, lowshelf_G, 1/sqrt(2), 1/sqrt(2), 'III');
    ls     = freqz(b, a, N/2);
else
    ls = ones(1,N/2)';
end

% HIGHSHELVING
if highshelf_apply
    if verbose
        disp('high shelving....');
    end
    [b, a] = AKhighshelve2(highshelf_f, fs, highshelf_G, 1/sqrt(2), 1/sqrt(2), 'III');
    hs     = freqz(b, a, N/2);
else
    hs = ones(1,N/2)';
end

% APPLY NOISE FLOOR
if length(noise_ir)~= 1 %#ok<*NODEF>
    % check input format
    if ischar(noise_ir)
        noise_apply = noise_ir;
        noise_ir = audioread(noise_ir);         % audioread if string
    else
        noise_apply = 1;
    end
    if size(noise_ir,2) > size(noise_ir,1)
        noise_ir = noise_ir';                 % transpose
    end
    if size(noise_ir,2)~=1
        noise_ir = mean(noise_ir, 2);         % average channels
    end
    if length(noise_ir) > N                   % adjust length
        noise_ir = noise_ir(1:N);
    elseif length(noise_ir) < N
        noise_ir(end:N) = 0;
    end
    % fft
    noise_ir = fft(noise_ir);
    noise_ir = noise_ir(1:N/2);
    % smoothing (only abslolute left after smoothing)
    if noise_smoothing
        noise_ir = AKfractOctSmooth(noise_ir, 'welti', fs, noise_smoothing);
    end
    % dynamic restriction
    if noise_dynamic
        if any(highPass)
            f_min = round(highPass(1)/df)+1;
        else
            f_min = 1;
        end
        if any(lowPass)
            f_max = round(lowPass(1)/df)+1;
        else
            f_max = round((fs/2-df)/df)+1;
        end
        min_val = max(noise_ir(f_min:f_max)) * 10^(-abs(noise_dynamic)/20);
        noise_ir(f_min:f_max) = max(noise_ir(f_min:f_max), min_val);
        noise_ir(1:f_min)   = min_val;
        noise_ir(f_max:end) = min_val;
    end
    % normalization
    noise_ir = noise_ir/min(noise_ir);
    clearvars f_min f_max min_val
    % plot
    figure
    semilogx(0:df:(N/2-1)*df, 20*log10(noise_ir), 'k','LineWidth',1.2)
    xlabel 'f in Hz'; ylabel 'Magnitude in dB_{rel}'; title 'Processed noise floor magnitude response (is not saved)';
    set(gcf, 'color', [1 1 1])
else
    noise_ir = ones(1,N/2)';
    noise_apply = 0;
end

% APPLY LOUDSPEAKER FREQUENCY RESPONSE
if length(speaker_ir)~= 1
    % check input format
    if ischar(speaker_ir)
        speaker_apply = speaker_ir;
        speaker_ir = audioread(speaker_ir);     % audioread if string
    else
        speaker_apply = 1;
    end
    if size(speaker_ir,2) > size(speaker_ir,1)
        speaker_ir = speaker_ir';               % transpose
    end
    if size(speaker_ir,2)~=1
        speaker_ir = mean(speaker_ir, 2);       % average channels
    end
    if length(speaker_ir) > N                % adjust length
        speaker_ir = speaker_ir(1:N);
    elseif length(speaker_ir) < N
        speaker_ir(end:N) = 0;
    end
    % fft
    speaker_ir = fft(speaker_ir);
    speaker_ir = speaker_ir(1:N/2);
    % smoothing (only abslolute left after smoothing)
    if speaker_smoothing
        speaker_ir = AKfractOctSmooth(speaker_ir, 'welti', fs, speaker_smoothing);
    end
    % inversion
    speaker_ir = 1./speaker_ir;
    % dynamic restriction
    if speaker_dynamic
        if any(highPass)
            f_min = round(highPass(1)/df)+1;
        else
            f_min = 1;
        end
        if any(lowPass)
            f_max = round(lowPass(1)/df)+1;
        else
            f_max = round((fs/2-df)/df)+1;
        end
        max_val = min(speaker_ir(f_min:f_max)) * 10^(abs(speaker_dynamic)/20);
        speaker_ir = min(speaker_ir, max_val);

    end
    % normalization
    speaker_ir = speaker_ir/min(speaker_ir);
    clearvars f_min f_max min_val max_val
    % plot
    figure
    semilogx(0:df:(N/2-1)*df, 20*log10(speaker_ir), 'k','LineWidth',1.2)
    xlabel 'f in Hz'; ylabel 'Magnitude in dB_{rel}'; title 'Processed speaker magnitude response (is not saved)';
    set(gcf, 'color', [1 1 1])
else
    speaker_ir = ones(1,N/2)';
    speaker_apply = 0;
end

% MR = Magnitude Response:
MR = ls .* hs .* noise_ir .* speaker_ir .* h_sweep;

% BANDPASS
if any(highPass)
    HP = AKfilter(AKdirac(N), 'hp', highPass(1), 0, fs, highPass(2), 'butter');
    HP = fft(HP);
    for nn = 1:highPass(3)
        MR = MR .* HP(1:size(MR,1));
    end
end

if any(lowPass)
    LP = AKfilter(AKdirac(N), 'lp', lowPass(1), 0, fs, lowPass(2), 'butter');
    LP = fft(LP);
    for nn = 1:lowPass(3)
        MR = MR .* LP(1:size(MR,1));
    end
end

% get absolute values
sweep_abs = abs(MR);

%mirror halfvector up to fs-df and get absolute values
%sweep_abs = mirror_h2fs(abs(MR));
sweep_abs = [sweep_abs; 0; flipud(sweep_abs(2:end))];

%% -------------------------------------------- 5. caluclate phase response
if verbose
    disp('calculating group delay...');
end
%initialize group delay
tg=zeros(nyq, 1);

%groupdelay at nyquist frequency
if NFFT_double
    tg(nyq)=T/2-(t_gap/1000);
else
    tg(nyq)=T-(t_gap/1000);
end

%groupdelay at DC
tg(1)=0;
%groupdelay for first frequency bin
tg(2)=t_start/1000;

% FORMULA (11, p.40 )
sweep_power = sum(abs(sweep_abs(3:nyq).^2));
C = (tg(nyq)-tg(2))/sweep_power;

% FORMULA (10, p.40 )
for k=3:nyq
    tg(k)=tg(k-1)+C*abs(sweep_abs(k))^2;
end

% calculating phase from group delay
if verbose
    disp('calculating phase from group delay...');
end
sweep_ang=-cumsum(tg)*2*pi*(df);

%wrapping phase
if verbose
    fprintf('wrapping phase...\n');
end
sweep_ang=AKwrap(sweep_ang(1:nyq));

% check if phase is zero at nyquist frequency
if sweep_ang(nyq)~= 0
    if verbose
        fprintf('phase(nyq)=%2.2f, not ZERO, correction running... \n',sweep_ang(nyq));
    end
    
    %correcting new phase
    sweep_ang=AKcorrectPhase(sweep_ang(1:nyq),fs);
    
    %wrapping phase again??
    if(max(abs(sweep_ang))>pi)
        if verbose
            fprintf('wrapping phase...\n');
        end
        sweep_ang=AKwrap(sweep_ang);
    end
end

%mirroring phase up to fs-df
sweep_ang=[sweep_ang(1:nyq); flipud(-1*sweep_ang(2:nyq-1))];
sweep_ang(1)=pi;

%% -------------------------------------- 6. calculate sweep in time domain
%calculate the complex spectrum
SWEEP    = sweep_abs.*exp(1i*sweep_ang);
SWEEP(1) = abs(SWEEP(1));

%into time domain
sweep=ifft(SWEEP);

%normalize to avoid clipping when written to wav with 16bit (worst case; LSB = 2^-15)
sweep = sweep' / max(abs(sweep)) * (1-2^-15);

%cut the second half of sweep, if NFFT-Size was doubled while
%synthesizing
if NFFT_double
    sweep=sweep(1:N/2);
    N = N/2;
end

% Fading in and out, zeropadding to NFFT
N_gap  = round(t_gap/1000*fs);
N_in   = round(t_start/1000*fs);
N_crop = N-round(N_gap*0.9);
N_out  = round(min([.1*fs .1*N_gap]));       % N_gap*0.1; (Old code from A. Giese)

if N_crop < 0
    error('t_gap longer than sweep duration, consider a smaller value')
end

%fading in
fade_in       = sin(linspace(0, pi/2, N_in)).^2;
sweep(1:N_in) = sweep(1:N_in) .* fade_in;
%cropping
sweep = sweep(1:N_crop);
%fading out
fade_out      = cos(linspace(0, pi/2, N_out)).^2;
sweep(end-N_out+1:end) = sweep(end-N_out+1:end) .* fade_out;
%zeropadding
sweep(end+1:N) = 0;

sweep = sweep';

%% ----------------------------------------------- 7. envelope_manipulation
if env_manipulation
    if verbose
        disp('envelope manipulation...')
    end
    % get positions of maxima and minima of each sine period
    sweep_diff = sign(diff(sweep));
    sweep_diff = [sweep_diff(1); sweep_diff];
    
    sweep_diff = abs(sign(diff(sweep_diff)));
    sweep_diff = [sweep_diff; 0];
    idx = find(sweep_diff ~= 0);
    
    % calculate envelope through maxima
    sweep_env = abs(sweep(idx));
    sweep_env = interp1(idx, sweep_env, 1:N)';
    
    % discard values before get start and end
    sweep_env(1:floor(t_start/1000*fs)) = 1;
    sweep_env(ceil(N-t_gap/1000*fs):end) = 1;
    
    % calculate what frequency occurs at each sample
    f_tg = linspace(0, fs/2, length(tg));
    tg_N = round(tg*fs);
    frequencies = zeros(N,1);
    for m = 1:N-1
        frequencies(tg_N(m)+1:tg_N(m+1)+1) = f_tg(m+1);
    end
    
    % smooth
    at = 0;
    at  = 1-exp(-2.2/fs/at*1000);
    %rt = max(frequencies/2, 10);
    rt = frequencies/2;
    
    sweep_env_sm_prev = 0;
    % allocate memory for sweep_env_sm
    sweep_env_sm = zeros(size(sweep_env));
    % loop over samples
    for idx = 1:size(sweep_env, 1)
        if abs(sweep_env(idx)) >= sweep_env_sm_prev
            % attack
            sweep_env_sm(idx) = at*abs(sweep_env(idx)) + (1-at)*sweep_env_sm_prev;
        else
            % release
            rt_curr = 1-exp(-2.2/fs/rt(idx)*1000);
            sweep_env_sm(idx) = (1-rt_curr)*sweep_env_sm_prev;
        end
        
        sweep_env_sm_prev = sweep_env_sm(idx, :);
        
    end
    
    clearvars at rt idx sweep_env_sm_prev tg_N f_tg sweep_diff
    
    % apply enevelope to sweep
    sweep = sweep .* (1./sweep_env_sm);
end

%% -------------------------------------------- 8. assign output parameters
%normalize to avoid clipping when written to wav with 16bit (worst case; LSB = 2^-15)
sweep = sweep / max(abs(sweep)) * (1-2^-15);

f_group_delay = (0:df:fs/2)';

varargout = {sweep, tg, f_group_delay};

%% ---------------------------------------------- 9. save wav and write log
if wav_save
    if ~ischar(wav_save)
        error('wav_save must be a string')
    end
    
    [wav_path, wav_name, wav_ext] = fileparts(wav_save);
    if isempty(wav_ext)
        wav_ext = '.wav';
    end
    
    % if wav file is written, write log file and plot in any case!!!
    if ~ischar(log_save)
        log_save = fullfile(wav_path, wav_name);
        do_plot  = fullfile(wav_path, wav_name);
    end
    if verbose
        disp('wav and log file...')
    end
        
    % write wav file
    audiowrite(fullfile(wav_path, [wav_name wav_ext]), sweep, fs, 'BitsPerSample', wav_Nbit)
end

if log_save
    if ~ischar(log_save)
        error('log_save must be a string')
    end
    if verbose
        disp('saving log file...')
    end
    % parse wav_save string
    delimiter_pos = strfind(log_save, '.');
    if delimiter_pos
        log_save  = [log_save(1:delimiter_pos) 'txt'];
    else
        log_save = [log_save '.txt'];
    end
    
    % write log file
    fid = fopen(log_save, 'w');
    fprintf(fid, ['Sweep generated with sweep_synth.m on ' date '\n\n']);
    fprintf(fid, 'preset      = ''%-3s''\n\n', preset);
    fprintf(fid, 'fs          = %-5u \t sampling frequency\n', fs);
    fprintf(fid, 'NFFT        = %-5u \t sweep order\n', NFFT);
    fprintf(fid, 'NFFT_double = %-5u \t sweep order has been doubled for synthesis\n', NFFT_double);
    fprintf(fid, 't_start     = %-5u \t time gap in ms before sweep starts\n', t_start);
    fprintf(fid, 't_gap       = %-5u \t silence after sweep to wait for high frequency-decay\n', t_gap);
    fprintf(fid, '\n');
    if highPass
        fprintf(fid, 'Butterworth high pass (N=%d, f_c=%d Hz) applied %d times to restrict sweep bandwidth\n', highPass(2), highPass(1), highPass(3));
    else
        fprintf(fid, 'No high pass applied\n');
    end
    if lowPass
        fprintf(fid, 'Butterworth low pass (N=%d, f_c=%d Hz) applied %d times to restrict sweep bandwidth\n', lowPass(2), lowPass(1), lowPass(3));
    else
        fprintf(fid, 'No low pass applied\n');
    end
    fprintf(fid, '\n');
    fprintf(fid, 'lowshelf_apply = %-8u \t apply lowshelf for bass emphasis\n', lowshelf_apply);
    if lowshelf_apply
        fprintf(fid, 'lowshelf_f     = %-8.2f \t midgain frequency in Hz\n', lowshelf_f);
        fprintf(fid, 'lowshelf_G     = %-8.2f \t gain in dB\n', lowshelf_G);
    end
    fprintf(fid, '\n');
    fprintf(fid, 'highshelf_apply = %-8u \t apply highsehlf for high frequency emphasis\n', highshelf_apply);
    if highshelf_apply
        fprintf(fid, 'highshelf_f     = %-8.2f \t midgain frequency in dB\n', highshelf_f);
        fprintf(fid, 'highshelf_G     = %-8.2f \t gain in dB\n', highshelf_G);
    end
    fprintf(fid, '\n');
    if ischar(noise_apply)
        fprintf(fid, 'noise floor has been applied from wav-file: %s\n', noise_apply);
        fprintf(fid, 'noise_dynamic   = %-8.2f \t dynamic in dB\n', noise_dynamic);
        fprintf(fid, 'noise_smoothing = 1/%-8u \t fractional octave smoothing\n', noise_smoothing);
    elseif noise_apply == 1
        fprintf(fid, 'noise floor has been applied from matlab variable \n');
        fprintf(fid, 'noise_dynamic   = %-8.2f \t dynamic in dB\n', noise_dynamic);
        fprintf(fid, 'noise_smoothing = 1/%-8u \t fractional octave smoothing\n', noise_smoothing);
    else
        fprintf(fid, 'no noise floor has been applied \n');
    end
    fprintf(fid, '\n');
    if ischar(speaker_apply)
        fprintf(fid, 'speaker frequency response has been applied from wav-file: %s\n', speaker_apply);
        fprintf(fid, 'speaker_dynamic   = %-8.2f \t dynamic in dB\n', speaker_dynamic);
        fprintf(fid, 'speaker_smoothing = 1/%-8u \t fractional octave smoothing\n', speaker_smoothing);
    elseif speaker_apply == 1
        fprintf(fid, 'speaker frequency response has been applied from matlab variable \n');
        fprintf(fid, 'speaker_dynamic   = %-8.2f \t dynamic in dB\n', speaker_dynamic);
        fprintf(fid, 'speaker_smoothing = 1/%-8u \t fractional octave smoothing\n', speaker_smoothing);
    else
        fprintf(fid, 'no speaker frequency response has been applied \n');
    end
    fprintf(fid, '\n');
    if env_manipulation
        fprintf(fid, 'Envelope manipulation has been applied \n');
    else
        fprintf(fid, 'Envelope manipulation has been not applied \n');
    end
    fprintf(fid, '\n');
    if wav_save
        fprintf(fid, 'wav_save = %s \t initial location of wav-file\n', wav_save);
        fprintf(fid, 'wav_Nbit = %u \t Bit resolution of wav-file\n', wav_Nbit);
    else
        fprintf(fid, 'wav_save = %u \t no wav file was written\n', wav_save);
    end
    fprintf(fid, '\n');
    if ischar(do_plot)
        fprintf(fid, 'do_plot  = %s \t initial location of documentaion plot\n', wav_save);
    end
    fclose(fid);
        
end

%% -------------------------------------------------------------- 10. plots
%TIMESIGNAL AND SPECTRUM
if do_plot
    
    AKf(20,30)

    %plot the time-signal
    subplot(3,1,1);
    AKp(sweep, 't2d', 'xu', 's', 'lw', 1.2, 'fs', fs)
    title('Sweep time signal')
    
    %plot the spectrum
    subplot(3,1,2);

    sweep_desired = 20*log10(sweep_abs/max(sweep_abs));
    df            = fs/(length(sweep_desired));
    f_plot        = 0:df:df*(length(sweep_desired)-1);
    
    semilogx(f_plot,sweep_desired, 'color', [.7 .7 .7],'LineWidth',1.2);
    AKp(sweep, 'm2d', 'lw', 1.2, 'norm_d', 0, 'x', [max(1, fs/N) fs/2], 'dr', [-60 5], 'fs', fs)
    
    legend('desired', 'obtained', 'location', 'best')
    
    % plot highest possible reverberation time
    subplot(3,1,3)
    if NFFT_double
        semilogx(f_group_delay, T/2-tg, 'k', 'lineWidth', 1.2);
        axis([max(10, f_group_delay(2)) fs/2 0 ceil(T/2)])
    else
        semilogx(f_group_delay, T-tg, 'k', 'lineWidth', 1.2);
        axis([max(10, f_group_delay(2)) fs/2 0 ceil(T)])
    end
    
    set(gca, 'xTick', [10:10:100 200:100:1000 2000:1000:10000 20000 30000 40000], ...
             'xTickLabel', {10    '' '' '' '' '' '' '' '' ...
                            100   '' '' '' '' '' '' '' '' ...
                            '1k'  '' '' '' '' '' '' '' '' ...
                            '10k' '20k' '' '40k'})
    
    grid on; box on
    
    title({'Time until the sweep ends' '(make sure this is longer than the reverberation time)'});
    xlabel('frequency in Hz');
    ylabel('time in s');
    
    if ischar(do_plot)
        % parse do_plot string
        delimiter_pos = strfind(do_plot, '.');
        if delimiter_pos
            do_plot  = do_plot(1:delimiter_pos-1);
        end
        print('-dpdf', do_plot)
    end
    
end

if do_advanced_plot
    
    AKf(20,10)
   
    % plot group delay of non-harmonic distortion products
    if NFFT_double
        AKsweepDistortion(tg, 2*N, fs);
    else
        AKsweepDistortion(tg, N, fs);
    end
    
    if ischar(do_plot)
        % parse do_plot string
        delimiter_pos = strfind(do_plot, '.');
        if delimiter_pos
            do_plot  = do_plot(1:delimiter_pos-1);
        end
        print('-dpdf', [do_plot ' distortionProducts'])
    end

end
