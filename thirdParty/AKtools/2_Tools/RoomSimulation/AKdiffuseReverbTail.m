% [reverb_tail, T_interp, f_interp] = AKdiffuseReverbTail(T, f, varargin)
%
% Models an n-channel diffuse reverbaration tail with frequency dependend
% decay based on reverberation time, binaural coherence, and diffuse field
% transfer function following [1], and [2].
%
% For examples see
% AKroomSimulationDemo.m, and
% AKdiffuseReverbTailDemo.m
%
%
% I N P U T:
% T             - Reverberation time in s
%                 (single value or vector with one row)
% f             - Frequencies in Hz of reverberation times given in T.
%                 f MUST increase monoton!
%
% OPTIONAL:
% V             - Volume of the room that you want to model (default = 200 m^3)
% c             - Speed of sound (default = 343 m/s)
% fs            - Sampling frequency (default = 44100 Hz)
% rev_channel   - Number of channles (default = 2)
% extrap_mode   - Method for extrapolating the reverbereation times below,
%                 inside and above tha range given by f.
%                 E.g. {'nearest', 'spline', 'linear'} which is the default
%                 takes the nearest neighbor below f(1), performs
%                 spline interpolation between f(1) and f(end), and
%                 performs spline interpolation above f(end). Possible
%                 values 'nearest', 'spline', 'linear'. Interpolation is
%                 done on a log(f) scale.
% extrap_lim    - Upper and lower limit for T (default = [8 0.01]). If these
%                 limits are exceeded by the interpolation, T will be
%                 clipped accordingly
% decay_mode    - Method for applying the exponential decay on the
%                 reverbration tail
%
%                 'non_stat' (default) uses non-stationary combination with
%                 a frequency resolution of T given by fs/N
%                 -> high resolution in frequency and time, long processing
%                 time
%
%                 'filt_bank' using a classical 3rd octave filter bank
%                 approach to apply one T for each filter band. Uses
%                 brickwall fft filter-bank from AKfilter.m
%                 -> High time resolution, limited frequency resolution,
%                 medium processing time
%
%                 'fft' applies frequency dependend decay per non
%                 overlappingfft block of length N
%                 -> Adjustable trade of between time and frequency
%                    resolution, shortest processing time
%
%                 'fft_ola' uses overlap and add for a smoother decay
%                 across fft blocks
%                 -> Adjustable trade of between time and frequency
%                    resolution, second shortest processing time
%
% decay_range   - Specifies the desired decay in dB up to which the tail is
%                 calculated (default = 90 dB). E.g. for max(T)=1 and
%                 decay_range=90, the reverb_tail will be 1.5 s long.
% bin_coherence - Applies binaural correlation based on spherical head
%                 model as introduced in [2] (default = true)
% dtf           - Applies a diffuse fied transfer function. Must be a
%                 single or, two channel impulse response (one, or two
%                 columns). Pass false to not apply a dtf (default = false)
% N             - Lenght of non-stationary filter for applying the decay if
%                 using 'non_stat', and length of fft block size if using
%                 'fft' (default = 2^10)
% N_window      - Length of window to fade out the reverb_tail
%                 (default = 2^5)
% do_plot       - Reverberation tail, interpolated T, and decay vs time and
%                 frequency (default = true)
%
%
% O U T P U T:
% reverb_tail   - n-channel reverberation tail
% T_interp      - inter/extrapolated reverberation time [seconds]
% f_interp      - frequency vector corresponding to values in T_interp [hertz]
%
%
% [1] James A Moorer (1979): "About this reverberation business." Computer
%     Music Journal, 3(2):13-28.
% [2] Christian Borss and Rainer Martin (2009): "An Improved Parametric
%     Model for Perception-Based Design of Virtual Acoustics." In: Proc.
%     of the 35th International AES Conference: Audio for Games. London
%
%
% PROCESSING:
% 1.  parse input
% 2.  generate noise and apply decay
% 3.  apply diffuse field HRTF
% 4.  binaural coherence
% 5.  normalize and fade out
% 6.  plot
%
% Audio Communication Group, TU Berlin, 11/2014
% fabian.brinkmann@tu-berlin.de; alexander.lindau@tu-berlin.de

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
function [reverb_tail, T_interp, f_interp] = AKdiffuseReverbTail(T, f, varargin)

%% ------------------------------------------ 1a. define default parameters
inputnames = {...
    'V', 'c', 'fs', 'rev_channel', ...
    'extrap_mode', 'extrap_lim', ...
    'decay_mode', 'decay_range', ...
    'bin_coherence', 'dtf'...
    'N', 'N_window', ...
    'do_plot'};

% Define default values
def.V             = 200;
def.c             = 343;
def.fs            = 44100;
def.rev_channel   = 2;
def.extrap_mode   = {'nearest' 'spline' 'linear'};
def.extrap_lim    = [8 .01];
def.decay_mode    = 'non_stat';
def.decay_range   = 90;
def.bin_coherence = true;
dtf               = false;
def.N             = 2^10;
def.N_window      = 2^5;
def.do_plot       = 1;

% --------------------------------------------------------- 1b. parse input
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

if bin_coherence && rev_channel~=2
    warning('AKdiffuse_reverb_tail:input', '''bin_coherence'' was set to false because ''rev_channel'' is not 2')
    bin_coherence = false;
end

% check format of common transfer function
% if ~islogical(dtf)
%     if size(dtf, 2) ~=2
%         dtf = [dtf(:,1) dtf(:,1)];
%     end
% end


%% --------------------------------------------------------- 2. apply decay
switch decay_mode
    case 'non_stat'
        % --- interpolate reverberation times ---
        % in case of non-stationary filtering for generation of the reverb, we want
        % to avoid having a bin at fs/2 in our fft. This bin had to be real and
        % would commonly be set to 0 (-inf dB) and thus we would have a dip in the
        % spectrum...
        if ~mod(N, 2)
            N = N+1;
        end
        
        % frequency vector
        f_interp = (0:fs/N:fs/2)';
        
        % get interpolated reverberation time
        T_interp = AKinterpolate_T(T, f, f_interp, fs, N, extrap_mode, extrap_lim);
        
        % --- generate the noise ---
        % get length of binaural tail in samples
        T_max = ceil(decay_range/60 * max(T_interp) * fs);
        
        % time axis
        t = (0:T_max-1)/fs;
        
        % get gaussian white noise
        reverb_noise = randn(T_max,rev_channel);
        
        % initialize reverberation tail
        reverb_tail   = zeros(size(reverb_noise));
        noise_channel = size(reverb_noise,2);
        
        % --- get decay ---
        % delta from decay formula in Heinrich Kutruff (2009): Room acoustics
        % 5.ed. Oxon: Spon Press, chpt. 5
        A     = 0.161*V./T_interp;
        delta = c*A/8/V;
        
        is_even = 1-mod(N,2);
        
        % save 100 impulse responses for plotting
        n_save    = floor(T_max/100);
        id_save   = 1;
        decay_course = zeros(N, 100);
        
        h = waitbar(0,'Waiting sucks - 0% done');
        
        % --- apply decay to noise ---
        for n=1:T_max
            % generate impulse response from decay values across frequency
            decay    = sqrt(exp(-2*delta.*t(n)));
            decay    = AKsingle2bothSidedSpectrum(decay, is_even);
            decay    = ifft(decay, 'symmetric');
            decay    = circshift(decay, round(N/2));
            decay    = phase_manipulation(decay, fs, 'min', 1, 0);
            
            % non-stationary combination
            nb = min([N n]);     % upper limit of hrir samples used
            nm = max([1 n-N+1]); % lower limit of input signal samples
            for ch = 1:noise_channel
                reverb_tail(n,ch) = fliplr(decay(1:nb)') * reverb_noise(nm:n,ch);
            end
            
            if ~mod(n-1, n_save)
                decay_course(:, id_save) = decay;
                waitbar((id_save-1)/100,h, ['Waiting sucks - ' num2str(id_save-1) '% done.'])
                id_save               = id_save + 1;
            end
        end
        
        % save t_course for plotting the course of the decay - similar to
        % waterfall plot
        t_course = t(1:n_save:end);
        
        %close waitbar
        close(h)
    case 'filt_bank'
        % --- interpolate reverberation times ---
        
        % frequency vector for interpolating reverberation time
        f_interp = (0:22050)';
        
        % get interpolated reverberation time
        T_interp = AKinterpolate_T(T, f, f_interp, fs, 44100, extrap_mode, extrap_lim);
        
        % --- generate the noise ---
        % get length of binaural tail in samples
        T_max = ceil(decay_range/60 * max(T_interp) * fs);
        
        % time axis
        t = (0:T_max-1)/fs;
        
        % get gaussian white noise
        reverb_tail = randn(T_max,rev_channel);
        % filter it
        [reverb_tail, fc] = AKfilter(reverb_tail, 'FFTfracOct', [0 fs/2], 0, 44100, repmat(3, [1 rev_channel]));
        
        % T at center frequencies of filtered noise
        T_interp = T_interp(round(fc(:,1))+1);
        f_interp = f_interp(round(fc(:,1))+1);
        
        % --- get decay ---
        % delta from decay formula in Heinrich Kutruff (2009): Room acoustics
        % 5.ed. Oxon: Spon Press, chpt. 5
        A     = 0.161*V./T_interp;
        delta = c*A/8/V;
        delta = repmat(delta', [T_max 1]);
        
        decay    = sqrt(exp(-2*delta.*repmat(t', [1 size(delta, 2)])));
        
        % apply decay
        for n = 1:rev_channel
            reverb_tail(:,n,:) = squeeze(reverb_tail(:,n,:)) .* decay;
        end
        
        % sum filterbands
        reverb_tail = sum(reverb_tail, 3);
        
    case 'fft'
        % --- interpolate reverberation times ---
        % in case of non-stationary filtering for generation of the reverb, we want
        % to avoid having a bin at fs/2 in our fft. This bin had to be real and
        % would commonly be set to 0 (-inf dB) and thus we would have a dip in the
        % spectrum...
        if ~mod(N, 2)
            N = N+1;
        end
        
        % frequency vector
        f_interp = (0:fs/N:fs/2)';
        
        % get interpolated reverberation time
        T_interp = AKinterpolate_T(T, f, f_interp, fs, N, extrap_mode, extrap_lim);
        
        % --- generate the noise ---
        % get length of binaural tail in samples
        T_max = ceil(decay_range/60 * max(T_interp) * fs);
        % adjust to be multiple of block length
        T_max = T_max + N-rem(T_max, N);
        % get number of blocks
        M = T_max/N;
        
        % get gaussian white noise
        reverb_tail = fft(randn(N, M, rev_channel));
        
        % --- get decay ---
        % delta from decay formula in Heinrich Kutruff (2009): Room acoustics
        % 5.ed. Oxon: Spon Press, chpt. 5
        A     = 0.161*V./T_interp;
        delta = c*A/8/V;
        
        % time axis
        t = (0:M-1)*N/fs;
        
        decay = sqrt(exp(-2*repmat(delta, [1 numel(t)]).*repmat(t, [numel(delta) 1])));
        decay = AKsingle2bothSidedSpectrum(decay, N);
        
        % apply decay
        for n = 1:rev_channel
            reverb_tail(:,:,n) = reverb_tail(:,:,n) .* decay;
        end
        
        % reshape
        reverb_tail = ifft(reverb_tail, 'symmetric');
        reverb_tail = reshape(reverb_tail, [T_max, rev_channel]);
        
        
    case 'fft_ola'
        % --- interpolate reverberation times ---
        % in case of non-stationary filtering for generation of the reverb, we want
        % to avoid having a bin at fs/2 in our fft. This bin had to be real and
        % would commonly be set to 0 (-inf dB) and thus we would have a dip in the
        % spectrum...
        if mod(N, 2)
            N = N+1;
        end
        
        % frequency vector
        f_interp = (0:fs/N:fs/2)';
        
        % get interpolated reverberation time
        T_interp = AKinterpolate_T(T, f, f_interp, fs, N, extrap_mode, extrap_lim);
        
        % --- generate the noise ---
        % get length of binaural tail in samples
        T_max = ceil(decay_range/60 * max(T_interp) * fs);
        % adjust to be multiple of block length
        T_max = T_max + N-rem(T_max, N);
        % get number of (non overlapping) blocks
        M = T_max/N;
        
        % get gaussian white noise
        reverb_noise = randn(N*M, rev_channel);
        
        % cut reverb noise into 2*M-1 blocks of size N with hop-size N/2
        reverb_tail               = zeros(N, 2*M-1, rev_channel);
        reverb_tail(:, 1:2:end,:) = reshape(reverb_noise,                   [N, M,   rev_channel]);
        reverb_tail(:, 2:2:end,:) = reshape(reverb_noise(N/2+1:end-N/2, :), [N, M-1, rev_channel]);
        
        reverb_tail = fft(reverb_tail);
        
        % --- get decay ---
        % delta from decay formula in Heinrich Kutruff (2009): Room acoustics
        % 5.ed. Oxon: Spon Press, chpt. 5
        A     = 0.161*V./T_interp;
        delta = c*A/8/V;
        
        % time axis for 2*M-1 blocks of size N with hop size N/2
        t = (0:.5:M-1)*N/fs;
        
        % actual decay
        decay = sqrt(exp(-2*repmat(delta, [1 numel(t)]).*repmat(t, [numel(delta) 1])));
        decay = AKsingle2bothSidedSpectrum(decay, N);
        
        % apply decay
        for n = 1:rev_channel
            reverb_tail(:,:,n) = reverb_tail(:,:,n) .* decay;
        end
        
        % --- overlap and add ---
        reverb_tail = ifft(reverb_tail, 'symmetric');
        
        % generate and applay window
        win              = sin(linspace(0, pi, N)').^2;
        win              = repmat(win, [1 2*M-1, rev_channel]);
        win(1:N/2,1,:)   = 1;
        win(N/2+1:end,end,:) = 1;
        
        reverb_tail = reverb_tail .* win;
        
        % overlap and add
        tmp_a = reshape(reverb_tail(:,1:2:end,:), [N*M,     rev_channel]);
        tmp_b = reshape(reverb_tail(:,2:2:end,:), [N*(M-1), rev_channel]);
        
        reverb_tail                    = tmp_a;
        reverb_tail(N/2+1:end-N/2,:,:) = reverb_tail(N/2+1:end-N/2,:,:) + tmp_b;
        
end

%% ------------------------------------------- 3. applay diffuse field HRTF
if any(dtf)
    reverb_tail = fftfilt(dtf, reverb_tail);
end

%% --------------------------------------------------------- 4. decorrelate
if bin_coherence
    reverb_tail = AKbinauralCoherence(reverb_tail, fs, N);
end

%% ---------------------------------------------- 5. normalize and fade out
reverb_tail = reverb_tail / max(abs(reverb_tail(:)));
if N_window
    reverb_tail(end-N_window+1:end, :) = reverb_tail(end-N_window+1:end, :) .* repmat(cos(linspace(0, pi/2, N_window)'),[1 size(reverb_tail,2)]);
end

%% ---------------------------------------------------------------- 6. plot
if do_plot
    
    if ~exist('decay_course', 'var')
        % make a spektogramm of the decay
        decay_course = reverb_tail(:,1);
        
        N = floor(numel(decay_course)/2^10);
        
        decay_course = reshape(decay_course(1:N*2^10), [2^10 N]);
        decay_course = ifft(abs(fft(decay_course)) ./ repmat(abs(fft(decay_course(:,1))), [1 N]), 'symmetric');
        
        t_course = (0:N) * 2^10/fs;
    end
    
    AKf(30,20)
    subplot(2,2,1)
        AKp(reverb_tail, 'et2d', 'c', 'cyc', 'dr', [-decay_range 0])
        title('ETC of reverb tail')
    subplot(2,2,2)
        AKp(reverb_tail, 'm2d', 'c', 'cyc')
        title('Magnitude spectrum of reverb tail')
    subplot(2,2,3)
        plot(f_interp,  T_interp, 'k')
        set(gca, 'xTick', f, 'xTickLabel', f/1000, 'xScale', 'log')
        xlabel('f in kHz'); ylabel('T in s')
        xlim([f_interp(1) 20000])
        title('Reverberation time')
    subplot(2,2,4)
        AKp(decay_course, 's3d', 'y', t_course, 'dr', [-decay_range 0])
        ylabel('t in s')
        title('Energy decay vs. frequency and time')
end

end



%% T_interp = AKinterpolate_T(T, f, f_interp, fs, N, extrap_mode, extrap_lim)

% T_interp = AKinterpolate_T(T, f, f_interp, fs, N, extrap_mode, extrap_lim)
% interpolates values at sparse frequencies to a frequency vector with
% fixed distances between neighboring frequencies
%
% INPUT:
% T           - input data [N x 1]
% f           - frequency bins corresponding to values of T
% f_interp    - frequencies to interpolate to
% fs          - sampling frequency in Hz, satisfying f<=fs/2
% N           - length of spectrum
% extrap_mode - Method for extrapolating the reverbereation times below,
%               inside and above tha range given by f. 
%               E.g. {'nearest', 'spline', 'spline'} which is the default
%               takes the nearest neighbor below f(1), performs
%               spline interpolation between f(1) and f(end), and
%               performs spline interpolation above f(end). Possible
%               values 'nearest', 'spline', 'linear'. Interpolation is
%               done on a log(f) scale.
% extrap_lim  - Upper and lower limit for T (default = [20 0.005]). If
%               these limits are exceeded by the interpolation, T will be
%               clipped accordingly

function T_interp = AKinterpolate_T(T, f, f_interp, fs, N, extrap_mode, extrap_lim)

% get bin of lowest and highest specified T values
f_lim = [ceil(f(1)/fs*N+1) floor(f(end)/fs*N+1)];

if numel(T) > 1
    % interpolate reverberation times
    T_interp                      = zeros(size(f_interp));
    T_interp(f_lim(1):f_lim(2)) = interp1(log(f), T, log(f_interp(f_lim(1):f_lim(2))), extrap_mode{2});
    
    % extrapolate reverberation times
    T_interp(1:f_lim(1)-1) = interp1(log(f), T, log(f_interp(1:f_lim(1)-1)), extrap_mode{1}, 'extrap');
    T_interp(f_lim(2)+1:end) = interp1(log(f), T, log(f_interp(f_lim(2)+1:end)), extrap_mode{3}, 'extrap');
    
    % clip reverberation times
    T_interp = min(T_interp, extrap_lim(1));
    T_interp = max(T_interp, extrap_lim(2));
    
    % 0 Hz might be critical when using interpolation on a log(f)
    % x_axis
    T_interp(1) = T_interp(2);
else
    T_interp = ones(size(f_interp))*T;
end

end