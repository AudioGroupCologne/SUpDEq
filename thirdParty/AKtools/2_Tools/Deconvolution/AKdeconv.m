% h = AKdeconv(y, x, fs, varargin)
%
% spectral deconvolution of data y with excitation signal x to obtain the
% impulse response h:
% H = Y/X (eq. 1),
% where H, Y, X are complex spectra, and h = ifft(H) (eq. 2).
% The division is in eq. 1 is done for each bin.
%
% h = AKdeconv(y, x, fs)
% will run the deconvolution using the default parameters.
% Additional data is passed in 'attribute' value pairs.
%
% For more examples see: AKdeconvolutionDemo.m
%
%
% I N P U T (default values):
% y                 - input time signal (samples x channels). Mono or
%                     multi-channel. This is the recorded sweep signal
% x                 - excitation signal (samples x channels). Mono or same
%                     number of channels as y
% fs                - sampling frequency in Hz
% deconv_type (lin) - 'cyc' for cyclic or 'lin' for linear deconvolution
%                     In linear deconvolution non harmonic distortion
%                     products are cut from the impulse response. This is
%                     done by doubling the length of x, and y by zero
%                     padding.
% x_inv_dyn (60)    - Inversion dynamic applied to x in dB: 1/X (eq. 1) is
%                     clipped if the desired dynamic range is exceeded.
%                     This might be usefull in the case of band limited
%                     excitation signals. Otherwise high gains in 1/X close
%                     to 0 Hz or fs/2 Hz can boost the noise in Y and
%                     influence the signal to noise ratio in H
% sys_comp (0)      - data for measurement system compensation, e.g.
%                     speaker or electroacoustical setup. sys_comp is
%                     inverted and applyed to x.
%                     format: matlab vector (samples x channels), .mat,
%                     .wav or .spk. Number of channels same as x or
%                     single channel.
% sys_inv_dyn (60)  - Inversion dynamic applied to sys_comp in dB
% sys_inv_phase (0) - phase that is applied to sys_inv after inversion:
%                     0 (no phase manipulation); 'min', 'lin'; 'zero'
% mic_comp (0)      - data for microphone compensation (see sys_comp)
% mic_inv_dyn (60)  - see sys_inv_dyn
% mic_inv_phase (0) - see sys_inv_phase
% highPass (false)  - restricts band width of 1/X (eq. 1) with a
%                     Butterworth high pass. See x_inv_dyn for reasons.
%                     Three element vector:
%                     [cut-off frequency in Hz, Order, M]
%                     The high pass is applied M times. Pass false, to
%                     disable.
% lowPass (false)   - restricts band width of 1/X (eq. 1) with a
%                     Butterworth low pass (c.f. highPass)
% do circshift (0)  - number of samples y is circshiftet to remove or add delay
% h_trunc (0)       - length in samples for shortening after deconvolution
% h_fade_in (0)     - length in samples for fade out (cos^2)
% h_fade_out (0)    - length in samples for fade in (cos^2)
% plot_deconv (0)   - plots results of deconvolution
% plot_filter (0)   - plots high-pass and low-pass filters
% plot_inv(0)       - plots results of inverting excitation signal, 
%                     sys_comp (if available) and mic_comp (if available)
%
% O U T P U T:
% h              - deconvolved time signal (samples x channels)
%
% PROCESSING
% 0: Define defualt parameters, data handling
% 1: Adjust length by zeropadding
% 2: Calculate spectra and inverse spectra
% 3: Filter inverse spectra with high- and lowpass
% 4: Spectral division and ifft
% 5: circshift and cut to original length if lin deconv
% 6: Truncate and window deconvolved ir
% 7: Plots
%
% issues: multichannel excitation signals have not been tested
%
% (c) F. Brinkmann, fabian.brinkmann@tu-berlin.de, 
% Audio Communication Group, TU Berlin
% 05/2012 v0.2 initial development
% 11/2012 V0.3 changed inversion dynamic to maintain phase information,
%              added plot of inverted excitation signal

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
function h = AKdeconv(y, x, fs, varargin)

%% ------------------------------------------- 0. define default parameters
inputnames = {...
'deconv_type',...
'x_inv_dyn',...
'sys_comp', 'sys_inv_dyn', 'sys_inv_phase',...
'mic_comp', 'mic_inv_dyn', 'mic_inv_phase',...
'highPass', ...
'lowPass', ...
'do_circshift', ...
'h_trunc', 'h_fade_in', 'h_fade_out',...
'plot_deconv', 'plot_filter', 'plot_inv'};

% Define default values
def.deconv_type    = 'lin';

def.x_inv_dyn      = 60;

def.sys_comp       = 0;
def.sys_inv_dyn    = 60;
def.sys_inv_phase  = 0;

def.mic_comp       = 0;
def.mic_inv_dyn    = 60;
def.mic_inv_phase  = 0;

def.highPass       = false;
def.lowPass        = false;

def.do_circshift   = 0;

def.h_trunc        = 0;
def.h_fade_in      = 0;
def.h_fade_out     = 0;

def.plot_deconv    = 0;
def.plot_filter    = 0;
def.plot_inv       = 0;

%% --------------------------------------------------------- 0. parse input
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


% ------------------------- congrats! you reached it to the processing part 
% Length and number of channels
N = max([size(y, 1) size(x, 1)]);
C = size(y, 2);

%% ---------------------------------------- 1. Adjust length by zeropadding
if strcmpi(deconv_type, 'lin')
    n = ceil(log2(N)) + 1;
    N = 2^n;
elseif ~strcmpi(deconv_type, 'cyc')
    error('Deconvolution deconv_type not defined, try cyc or lin')
end

% zero pad
y(end+1:N, :) = 0;
x(end+1:N, :) = 0;

%% ------------------------------- 2. calculate spectra and inverse spectra
% FFT und inverte reference H
Y = fft(y, [], 1);
X = fft(x, [], 1);

% ------ Inverse excitation signal ------
X_inv = 1 ./ X;
% apply inversion dynamic
if x_inv_dyn
    tmp_min = min(min(abs(X_inv)));
    X_inv(abs(X_inv) > tmp_min*10^(abs(x_inv_dyn)/20)) = ...                % select bins that exceed inversion dynamic,
        tmp_min*10^(abs(x_inv_dyn)/20) .* ...                               % substitute magnitude and
        exp(1j*angle(X_inv(abs(X_inv) > tmp_min*10^(abs(x_inv_dyn)/20))));  % leave phase untouched
end
% data format
if size(X_inv, 2) ~= C
    X_inv = repmat(X_inv, [1 C]);
end

X_inv_all = X_inv;

% ----- inverse microfone transfer function -----
if max(mic_comp)
    % load, truncate or zerp pad and get inverse
    if ischar(mic_comp)
        mic_comp = data_read(mic_comp);
    end
    if length(mic_comp) > N
        mic_comp = mic_comp(1:N, :);
    else
        mic_comp(end+1:N, :) = 0;
    end
    Mic_Comp_Inv = 1./fft(mic_comp);
    % restric inversion dynamic
    if mic_inv_dyn
        tmp_min = min(min(abs(Mic_Comp_Inv)));
        Mic_Comp_Inv(abs(Mic_Comp_Inv) > tmp_min*10^(abs(mic_inv_dyn)/20)) = ...                % select bins that exceed inversion dynamic,
            tmp_min*10^(abs(mic_inv_dyn)/20) .* ...                                             % substitute magnitude and
            exp(1j*angle(Mic_Comp_Inv(abs(Mic_Comp_Inv) > tmp_min*10^(abs(mic_inv_dyn)/20))));  % leave phase untouched
    end
    % manipulate phase information
    if mic_inv_phase
        Mic_Comp_Inv = real(ifft(Mic_Comp_Inv));
%         Mic_Comp_Inv = make_phase(Mic_Comp_Inv, mic_inv_phase);
        Mic_Comp_Inv = AKphaseManipulation(Mic_Comp_Inv, fs, mic_inv_phase, 1, 0);
        Mic_Comp_Inv = fft(Mic_Comp_Inv);
    end
    % data format
    if size(Mic_Comp_Inv, 2) ~= C
        Mic_Comp_Inv = repmat(Mic_Comp_Inv, [1 C]);
    end
    
    % combine inverse transfer functions
    X_inv_all = X_inv_all .* Mic_Comp_Inv;
end

% ----- inverse measurement system transfer function -----
if max(sys_comp)
    % load, truncate or zerp pad and get inverse
    if ischar(sys_comp)
        sys_comp = data_read(sys_comp);
    end
    if length(sys_comp) > N
        sys_comp = sys_comp(1:N, :);
    else
        sys_comp(end+1:N, :) = 0;
    end
    Sys_Comp_Inv = 1./fft(sys_comp);
    % restric inversion dynamic
    if sys_inv_dyn
        tmp_min = min(min(abs(Sys_Comp_Inv)));
        Sys_Comp_Inv(abs(Sys_Comp_Inv) > tmp_min*10^(abs(sys_inv_dyn)/20)) = ...                % select bins that exceed inversion dynamic,
            tmp_min*10^(abs(sys_inv_dyn)/20) .* ...                                             % substitute magnitude and
            exp(1j*angle(Sys_Comp_Inv(abs(Sys_Comp_Inv) > tmp_min*10^(abs(sys_inv_dyn)/20))));  % leave phase untouched
    end
    % manipulate phase information
    if sys_inv_phase
        Sys_Comp_Inv = real(ifft(Sys_Comp_Inv));
%         Sys_Comp_Inv = make_phase(Sys_Comp_Inv, sys_inv_phase);
        Sys_Comp_Inv = AKphaseManipulation(Sys_Comp_Inv, fs, sys_inv_phase, 1, 0);
        Sys_Comp_Inv = fft(Sys_Comp_Inv);
    end
    % data format
    if size(Sys_Comp_Inv, 2) ~= C
        Sys_Comp_Inv = repmat(Sys_Comp_Inv, [1 C]);
    end
    
    % combine inverse transfer functions
    X_inv_all = X_inv_all .* Sys_Comp_Inv;
end


%% ----------------------- 3. Filter inverse spectra with high- and lowpass
x_filt_HP      = zeros(N,C);
x_filt_HP(1,:) = 1;
x_filt_LP      = x_filt_HP;

if any(highPass)
    HP = AKfilter(AKdirac(N), 'hp', highPass(1), 0, fs, highPass(2), 'butter');
    for nn = 1:highPass(3)
        x_filt_HP = fftfilt(x_filt_HP, HP);
    end
else
    X_hp = 'No highpass applied';
end

if any(lowPass)
    LP = AKfilter(AKdirac(N), 'lp', lowPass(1), 0, fs, lowPass(2), 'butter');
    for nn = 1:lowPass(3)
        x_filt_LP = fftfilt(x_filt_LP, LP);
    end
else
    X_lp = 'No highpass applied';
end

% filter
X_filt = fft(x_filt_HP) .* fft(x_filt_LP);
X_inv_all_filt = X_inv_all .* X_filt;


%% ------------------------------------------ 4. Spectral division and ifft
% Divide spectra
H = Y .* X_inv_all_filt;

% IFFT
h = real(ifft(H));

%% ------------------ 5. circshift and cut to original length if lin deconv
if do_circshift
    h = circshift(h, [do_circshift, 0]);
end

% cut to original length
if strcmpi(deconv_type, 'lin')
    N = N/2;
    h = h(1:N, :);
end

%% ---------------------------------- 6. truncate and window deconvolved ir
% i.e. to discard artefacts at end of impulse response (distortions from
% deconvolution) or discard reflections
if plot_deconv
    h_pre_window = h;
end

% length for windowing
if h_trunc
    h = h(1:h_trunc, :);
end

window = ones(size(h));
if h_fade_in
    fade_tmp = sin(linspace(0, pi/2, h_fade_in)).^2;
    window(1:h_fade_in, :) = repmat(fade_tmp', [1 C]);
end
if h_fade_out
    fade_tmp = cos(linspace(0, pi/2, h_fade_out)).^2;
    window(end-h_fade_out+1:end, :) = repmat(fade_tmp', [1 C]);
end
clearvars fade_tmp

h = h .* window;

%% --------------------------------------------------------------- 7. plots
warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart') % p function complaines, though all data is real as shit
% ----- plot deconvolved signals -----
if plot_deconv
    for n = 1:C
        AKf(20,15)
        set(gcf, 'name', ['Deconvolved channel ' num2str(n)])
        % etc
        subplot(2,1,1)
            AKp(h_pre_window(:, n), 'et2d', 'norm_d', 0, 'c', .7, 'xu', 'n', 'fs', fs)
            AKp(h(:, n), 'et2d', 'norm_d', 0, 'xu', 'n', 'fs', fs)
            AKp(window(:, n), 'et2d', 'c', 'r', 'xu', 'n', 'fs', fs)
            title('Normalized ETC (black = final signal, red = window)')
            if h_trunc
                axis([1 h_trunc+500 -100 5])
            else
                axis([1 length(h) -100 5])
            end
        % spectrum
        subplot(2,1,2)
            AKp(h(:,n), 'm2d', 'fs', fs, 'x', [20 20000], 'c', .7, 'N', max(fs/20, numel(h(:,n))))
            AKp(h(:,n), 'm2d', 'fs', fs, 'x', [20 20000], 'frac', 3, 'N', max(fs/20, numel(h(:,n))))
            title('Magnitude response. Black = 3rd oct. sm.')
    end
end

% ----- plot filter -----
if plot_filter
    % plot hp
    if any(highPass)
        if any(highPass) && any(lowPass)
            AKf(20,10)
            subplot(1,2,1)
        else
            AKf(20,10)
        end
        % plot spectrum
        AKp(x_filt_HP, 'm2d', 'N', fs, 'lw', 1.2, 'x', [1 fs/2], 'fs', fs)
        title(['Butterworth high-pass, order:' num2str(highPass(2)) ', cut-off:' num2str(highPass(1)) ' Hz, (applied ' num2str(highPass(3)) ' times)'])
    end
    
    % plot lp
    if any(lowPass)
        if any(highPass) && any(lowPass)
            subplot(1,2,2)
        else
            AKf(20,10)
        end
        % plot spectrum
        AKp(x_filt_LP, 'm2d', 'N', fs, 'lw', 1.2, 'x', [1000 fs/2], 'fs', fs)
        title(['Butterworth low-pass, order:' num2str(lowPass(2)) ', cut-off:' num2str(lowPass(1)) ' Hz, (applied ' num2str(lowPass(3)) ' times)'])
    end
end

% ----- plot inverted microphone and system freq. response -----
if plot_inv
    
    AKf(15,8)
    NN = length(X_inv);
    f  = (0:NN-1)*fs/NN;
    semilogx(f, 20*log10(abs(X_inv.*X_filt)), 'k', 'linewidth', 1.2)
    set_xAxis_freq(gca)
    xlim([10 fs/2])
    xlabel('f/Hz'); ylabel('Magnitude/dB'); title('Inverted excitation signal freq. response')
    
    if exist('Mic_Comp_Inv', 'var')
        AKf(15,8)
        NN = length(Mic_Comp_Inv);
        f  = (0:NN-1)*fs/NN;
        semilogx(f, 20*log10(abs(Mic_Comp_Inv)), 'k', 'linewidth', 1.2)
        set_xAxis_freq(gca)
        xlim([10 fs/2])
        xlabel('f/Hz'); ylabel('Magnitude/dB'); title('Inverted Microphone freq. response')        
    end
    
    if exist('Sys_Comp_Inv', 'var')
        AKf(15,8)
        NN = length(Sys_Comp_Inv);
        f  = (0:NN-1)*fs/NN;
        semilogx(f, 20*log10(abs(Sys_Comp_Inv)), 'k', 'linewidth', 1.2)
        set_xAxis_freq(gca)
        xlim([10 fs/2])
        xlabel('f/Hz'); ylabel('Magnitude/dB'); title('Inverted System freq. response')        
    end
end
warning('on', 'MATLAB:plot:IgnoreImaginaryXYPart')

end

function set_xAxis_freq(gca)

 
grid on
set(gca, 'xTick', [1:1:10 20:10:100 200:100:1000 2000:1000:10000 20000:10000:100000 200000], ...
    'xTickLabel', {'1'    '' '' '' '' '' '' '' '' ...
                   '10'   '' '' '' '' '' '' '' '' ...
                   '100'  '' '' '' '' '' '' '' '' ...
                   '1k'   '' '' '' '' '' '' '' '' ...
                   '10k'  '' '' '' '' '' '' '' '' ...
                   '20k'  '' '' '' '' '' '' '' '' ...
                   '100k' '' '' '' '' '' '' '' '' ...
                   '200k'},...
     'GridLineStyle', '-')
end
