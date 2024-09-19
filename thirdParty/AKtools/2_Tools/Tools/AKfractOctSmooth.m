% X_s = AKfractOctSmooth(x,[operation,fs,frac,mode,band,b,do_plot]);
%
% Applies fractional octave band smoothing on multichannel spectrum
% (for WELTI delivered as absolute or complex single-sided spectrum)
% (for SCHAERER delivered as impulse responses)
%
% To smooth a noisy spectrum try this
% x = AKnoise(2^12);
% X = AKboth2singleSidedSpectrum(fft(x));
% Y = AKfractOctSmooth(X,'amp',44100,3,true);
%     % or
% Y = AKfractOctSmooth(x,'schaerer',44100,3,'amp','oct',[],true);
% y = ifft(AKsingle2bothSidedSpectrum(Y);
%
%
%  INPUT operation - 'amp' for pure amplitude smoothing (zero phase),
%                    'cmp' for seperate amplitude and phase smoothing,
%                    'eqv' for amplitude smoothing plus copying phase of
%                          input signal,
%                    default = 'amp'
%
%     output:
%        X_s       - smoothed single-sided spectrum with same length as the
%                    input spectrum
%     input:
%        x         - single or multi channel SINGLE-SIDED SPECTRA
%        fs        - sampling frequency, default = 44100
%        frac      - fraction of octave for smoothing,
%                    default = 3 (third octave smoothing)
%        do_plot   - [0 .. 1 .. channel number] show comparative plot,
%        [ fifth     given as true/false or number of channel which should
%          param ]   be shown, default = false
%
%                    ----------------------
%
%  INPUT operation - 'schaerer'
%                  - This algorithm implementation is quite expensive,
%                    especially for long input signals (> 4096 samples).
%                    Maybe apply signal truncation or windowing before!
%
%  O U T P U T:
%        X_s       - smoothed single-sided spectrum with HALF EVEN length
%                    as the input spectrum
%     input:
%        x         - single or multi channel IMPULSE RESPONSES
%        fs        - sampling frequency, default = 44100
%        frac      - fraction of octave for smoothing,
%                    default = 3 (third octave smoothing)
%        mode      - 'amp' for amplitude,
%                    'cpm' for complex,
%                    'eqcmp' for equivalent complex (complex with amplitude
%                            compensation),
%                    default = 'amp'
%        band      - 'oct' for fractional octave,
%                    'bark' for bark,
%                    'ERB' for equivalent rectangular bandwidth,
%                    default = 'oct'
%        b         - coefficient for smoothing window
%                    W = b - (1-b)*cos(2*pi*k/N), for k=[0,length(x)-1],
%                    i.e. .5 (Hanning), .54 (Hamming), 1 (rectangular)],
%                    default = .5
%        do_plot   - [0 / 1 / channel number] show comparative plot,
%                    given as true/false or number(s) of channel which 
%                    should be shown, default = false
%
% ----------------------------------------------------------------------- %
% Additional information to operational mode 'schaerer':
%
% Calculates the smoothed transfer function to an input [x]. Returns the
% smoothed transfer function X_s.
% The smoothed function is calculated by convolving the transfer function
% of H with a window function W. The width of the window increases with
% increasing frequency, [band] determines the frequency scale used
% (fraction of octave, bark or ERB). The convolution is either conducted
% with the magnitude transfer function or the complex transfer function. In
% the latter case, the window is also applied to the time domain,
% increasing in width with decreasing frequency. In the case of equivalent
% complex smoothing, X_s is a combination of the amplitude smoothed
% magnitude and the complex smoothed phase.
% [Hatziantoniou PD, Mourjopoulos JN (2000) Generalized Fractional-Octave
% Smoothing of Audio and Acoustic Responses]
% ----------------------------------------------------------------------- %
% v1    2009 Andreas Rotter, geeprombolo@gmail.com,
%            Audio Communicatin Group, TU Berlin,
%            using functions from Todd Welti and Andre Giese
% v2 01/2018 Hannes Helmholz, helmholz@campus.tu-berlin.de,
%            Audio Communicatin Group, TU Berlin,
%            rework and extention of Welti algorithm; improvement and
%            smoothing fraction adjustment of Schaerer algorithm
% v3 02/2018 fixed behaviour for real spectrum inputs (e.g. AKp)
% v4 02/2018 calculation time warning and better string handling
% v5 10/2019 Hannes Helmholz, hannes.helmholz@chalmers.se,
%            Applied Acoustics, Chalmers University of Technology,
%            improvement at frequency boundaries of Welti algorithm 
% ----------------------------------------------------------------------- %

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

function X_s = AKfractOctSmooth(x,operation,fs,frac,mode,band,b,do_plot)

% an arbitrary value over which a Matlab warning gets shown in case the
% approximated calculation time for the Schaerer algorithm excals the given
% value
CALC_TIME_WARNING_SEC = 30;

if nargin < 2 || isempty(operation)
    operation = 'amp';
end
if nargin < 3 || isempty(fs)
    fs = 44100;
end
if nargin < 4 || isempty(frac)
    frac = 3;
end
if strcmpi(operation,'schaerer')
    if nargin < 5 || isempty(mode)
        mode = 'amp';
    end
    if nargin < 6 || isempty(band)
        band = 'oct';
    end
    if nargin < 7 || isempty(b)
        b = .5;
    end
    if nargin < 8
        do_plot = 0;
    end
    
    switch lower(mode)
        case 'amp'
            fprintf('Schaerer -> amplitude smoothing');
        case 'cmp'
            fprintf('Schaerer -> complex smoothing');
        case 'eqcmp'
            fprintf('Schaerer -> equivalent complex smoothing');
        otherwise
            fprintf('\n');
            error('AKfractOctSmooth:mode','Illegal mode.');
    end
    fprintf(' -> ');
    switch lower(band)
        case 'oct'
            fprintf('bandwidth: 1/%g octave',frac);
        case 'bark'
            fprintf('bandwidth: critical band');
        case 'erb'
            fprintf('bandwidth: ERB');
        otherwise
            fprintf('\n');
            error('AKfractOctSmooth:band','Illegal band.');
    end
    fprintf(' -> ...');
else
    if nargin < 5
        do_plot = 0;
    else
        do_plot = mode;
    end
end

switch lower(operation)
    % WELTI is a leftover from an old version
    case {'amp','amplitude','abs','absolute','welti'}
        operation = 'amp';
    case {'cmp','complex','amp+pha','abs+pha','amp+ang','abs+ang'}
        operation = 'cmp';
    case {'eqv','equivalent','amp+cpy','abs+cpy'}
        operation = 'eqv';
    case 'schaerer'
        % nothing to do here
    otherwise
        error('AKfractOctSmooth:operation','Illegal operation.');
end

% skip smoothing
if ~frac
    if strcmpi(operation,'schaerer')
        X_s = AKboth2singleSidedSpectrum(fft(x));
        fprintf(' DONE.\n');
    else
        X_s = x;
    end
    return;
end

frac = 1/frac;

% ----------------------------------------------------------------------- %

if strcmpi(operation,'schaerer')
    time = SchaererApproxCalcTime(size(x,1),size(x,2));
    if time > CALC_TIME_WARNING_SEC
        fprintf('\n');
        warning('AKfractOctSmooth:schaerer',['The estimated '...
            'calculation time for %d impulse response(s) of %d samples '...
            'is about %d seconds.'],size(x,2),size(x,1),ceil(time));
        fprintf(' -> ...');
    end
    X_s = AKfractOctSmoothSchaerer(fft(x),fs,frac,mode,band,b);
    
    % ------------------------------------------------------------------- %
else
    % frequency support vector
    f = linspace(0,fs/2,size(x,1));
 
    if isreal(x) && strcmpi(operation,'amp')
        % direct smoothing, used e.g. by AKp
        X_s = AKfractOctSmoothWelti(x,frac,f(end),f);
        
    else
        % seperate magnitude and phase
        X_s_abs = AKfractOctSmoothWelti(abs(x).^2,frac,f(end),f);
        switch lower(operation)
            case 'amp'
                X_s_angle = zeros(size(x));
            case 'cmp'
                X_s_angle = AKfractOctSmoothWelti(angle(x),frac,f(end),f);
                % % One could do the smoothing on the unwrapped phase
                % % which leads to a different but not more meaningful IR.
                % % Panzer2004 �The Use of Continuous Phase for ..."
                % % suggests this should work somehow though.
                % % (H. Helmholz, 01/2018)
                % X_s_angle = AKfractOctSmoothWelti(unwrap(angle(x)),frac,f(end),f);
            case 'eqv'
                X_s_angle = angle(x);
        end
        X_s = sqrt(X_s_abs) .* exp(1i*X_s_angle);
    end
end

% ----------------------------------------------------------------------- %

% plotting
if do_plot
    channels = do_plot; % also the channel index which should be plotted
    if strcmpi(operation,'schaerer')
        operation_plot = [mode,' ',band];
        x_plot = x(:,channels);
    else
        operation_plot = operation;
        x_plot = ifft(AKsingle2bothSidedSpectrum(x(:,channels),true));
    end
    x_s_plot = ifft(AKsingle2bothSidedSpectrum(X_s(:,channels),true));
    
    AKf;
    set(gcf,'Name',sprintf('%s  smoothing (ch %s)',...
        operation_plot,AKlistIntegers(channels)),'NumberTitle','off');
    AKpMulti(x_plot,{'t2d','et2d';'m2d','p2d'},'c','k');
    AKpMulti(x_s_plot,{'t2d','et2d';'m2d','p2d'},'c','r','lw',2);
    AKtightenFigure;
end

if strcmpi(operation,'schaerer')
    fprintf(' DONE.\n');
end

end


% ----------------------------------------------------------------------- %
function H_smooth = AKfractOctSmoothWelti(H,oct,f_stop_i,f)

if nargin >= 4
    f_stop_i = round(f_stop_i/f(2))+1;
end

% cast H to double to avoid error (F. Brinkmann 7/2012)
H = cast(H,'double');

% fixed at 1 for now, doesn't work if not starting at 1
f_start_i = 1;
% number of samples for log warping, seems to work well
f_N = f_stop_i;
% this is a multiplicitive factor, not an even spacing
oct_spacing = 10^(log10(f_stop_i-f_start_i)/f_N);
% number of bins per fractional octave
oct_N = oct*log10(2)/log10(oct_spacing);

% logarithmic warping
f_log = logspace(log10(f_start_i),log10(f_stop_i),f_N);
H_log = interp1(1:size(H,1),H,f_log,'spline');
H_log = reshape(H_log,size(H));
% don't change to SPLINE(), because of interpolation per column (not row)

% make it even length for convenience
oct_N_even = round(oct_N/2)*2;
win = gausswin(oct_N_even*2);

% minimze smoothing error at boundaries by extension
H_ext = ones(oct_N_even,size(H_log,2),size(H_log,3));
% % VAR1: repeat value of first / last frequency bin
% H_begin = H_log(1,:,:) .* H_ext;
% H_end = H_log(end,:,:) .* H_ext;
% % VAR2: repeat value of average of half window length from first / last frequency bin
% H_begin = mean(H_log(1:oct_N_even,:,:),1) .* H_ext;
% H_end = mean(H_log(end+1-oct_N_even:end,:,:),1) .* H_ext;
% VAR3: repeat value of weighted average of half window length from first / last frequency bin
H_begin = sum(H_log(1:oct_N_even,:,:).*win(oct_N_even+1:end),1)./sum(win(oct_N_even+1:end));
H_begin = H_begin .* H_ext;
H_end = sum(H_log(end+1-oct_N_even:end,:,:).*win(1:oct_N_even),1)./sum(win(1:oct_N_even));
H_end = H_end .* H_ext;
% % VAR4: mirrorred values of half window length around first / last frequency bin
% H_begin = flip(H_log(1:size(H_ext,1),:,:),1);
% H_end = flip(H_log(end+1-size(H_ext,1):end,:,:),1);
H_log = [H_begin ; H_log ; H_end];

% convolve with smoothing window
H_smooth = zeros(size(H_log));
for ch = 1:size(H_log,ndims(H_log))
    H_smooth(:,ch) = fftfilt(win,H_log(:,ch),2*length(win))./sum(win);
end

% get rid of extra bins from adding lead and lag
H_smooth = H_smooth(size(H_ext,1)+length(win)/2 + (1:length(f_log)),:,:); 

% dewarping to linear
H_smooth = interp1(f_log,H_smooth,1:f_stop_i,'spline');
H_smooth = reshape(H_smooth,size(H));
% don't change to SPLINE(), because of interpolation per column (not row)

% transpose if needed to get same output as input (column/row)
if size(H,1)~=f_stop_i && size(H,2)~=1
    H_smooth = H_smooth';
end

end


% ----------------------------------------------------------------------- %
function H_smooth = AKfractOctSmoothSchaerer(H,fs,frac,mode,band,b)

% reshape matrix to also handle [filterlen,1,channels] inputs
H_orig_size = size(H);
H = squeeze(H);

[filterlen,channels] = size(H);
filterlen_2 = floor(filterlen/2);

if any(strcmpi(band,{'bark','erb'}))
    frac = 1/frac;
    b = frac;
end

% generate vector with bandwidths according to the selected scale
d_f = fs/filterlen;
switch lower(band)
    case 'oct'
        % % original Schaerer implementation
        % f_up = (10^(3/10)) ^ (.5*frac);
        % f_low = (10^(-3/10)) ^ (.5*frac);
        % bandwidth = (1:filterlen_2) * d_f*(f_up-f_low);
        
        % changes see below (H. Helmholz 01/2018)
        bandwidth = (1:filterlen_2) * d_f;
    case 'bark'
        bandwidth = 25+75*(1+1.4*((1:filterlen_2) * d_f/1000).^2) .^ .69;
    case 'erb'
        bandwidth = 24.7*(4.37*((1:filterlen_2) * d_f/1000) + 1);
end
% the sqrt(3*frac) correction factor was imcluded, to match the Welti
% algorithm smoothing behaviour (with b=.5 and band='oct')
% (H. Helmholz 01/2018)
m = ceil(.5 * bandwidth * sqrt(3*frac)/d_f);

% generate windowfunction
W = zeros(filterlen_2,filterlen);
for s = 1:filterlen_2
    M = m(s);
    W(s,1:M+1) = (b - (b-1)*cos((pi/M)*(0:M))) ./ (2*b*(M+1)-1);
    W(s,filterlen-M+1:filterlen) = ...
        (b - (b-1)*cos((pi/M)*((filterlen-M:filterlen-1) - filterlen))) ...
        ./ (2*b*(M+1)-1);
end

% seperate amplitde and phase
if any(strcmpi(mode,{'amp','eqcmp'}))
    H_2 = abs(H).^2;
    
    % generate circular convolution matrix
    H_circ_2 = zeros(filterlen,filterlen,channels);
    for s = 1:filterlen
        H_circ_2(s,:,:) = circshift(H_2,s-1,1);
    end
    
    % circular convolution
    H_smooth_circ_2 = zeros(filterlen_2,filterlen,channels);
    for ch = 1:channels
        H_smooth_circ_2(:,:,ch) = W * H_circ_2(:,:,ch);
    end
    
    % derive smoothed function
    H_smooth_2 = zeros(filterlen,channels);
    for s = 1:filterlen_2
        H_smooth_2(s,:) = squeeze(H_smooth_circ_2(m(s),s,:));
        if s > 1 % mirror spectrum
            H_smooth_2(filterlen-s+2,:) = H_smooth_2(s,:);
        end
    end
end
if any(strcmpi(mode,{'cmp','eqcmp'}))
    H_r = real(H);
    H_im = imag(H);
    
    % generate circular convolution matrix
    H_circ_r = zeros(filterlen,filterlen,channels);
    H_circ_im = H_circ_r;
    for s = 1:filterlen
        H_circ_r(s,:,:) = circshift(H_r,s-1,1);
        H_circ_im(s,:,:) = circshift(H_im,s-1,1);
    end
    
    % circular convolution
    H_smooth_circ_r = zeros(filterlen_2,filterlen,channels);
    H_smooth_circ_im = H_smooth_circ_r;
    for ch = 1:channels
        H_smooth_circ_r(:,:,ch) = W * H_circ_r(:,:,ch);
        H_smooth_circ_im(:,:,ch) = W * H_circ_im(:,:,ch);
    end
    
    % derive smoothed function
    H_smooth_r = zeros(filterlen,channels);
    H_smooth_im = H_smooth_r;
    for s = 1:filterlen_2
        H_smooth_r(s,:) = squeeze(H_smooth_circ_r(m(s),s,:));
        H_smooth_im(s,:) = squeeze(H_smooth_circ_im(m(s),s,:));
        if s > 1 % mirror spectrum
            H_smooth_r(filterlen-s+2,:) = H_smooth_r(s,:);
            H_smooth_im(filterlen-s+2,:) = H_smooth_im(s,:);
        end
    end
    H_smooth_c = H_smooth_r + 1i*H_smooth_im;
end

% combine amplitude and phase
switch lower(mode)
    case 'amp'
        H_smooth = sqrt(H_smooth_2) .* exp(1i*angle(H));
    case 'cmp'
        H_smooth = H_smooth_c;
    case 'eqcmp'
        H_smooth = sqrt(H_smooth_2) .* exp(1i*angle(H_smooth_c));
end

% restore input dimensions
H_smooth = reshape(H_smooth,H_orig_size);
H_smooth = AKboth2singleSidedSpectrum(H_smooth);

end


function time_sec = SchaererApproxCalcTime(ir_len,n_irs)
% time approximation based on these values on a 3.5 GHz 6-Core machine:
%
% len     = [1024,2048,4096,8192,16384];
% time    = [0.518675,0.92,3.2,16.205843,129.771364];
% weights = [.1,.25,.5,1,1];
%
% which leads to the following coefficient by curve fitting:

a = 1.31e-11;
b =    3.083;
c =    1.052;

time_sec =  n_irs * a * ir_len ^ b + c;

% This approximation will be not quite right for another machine of course.
% Also the scaling with an increasing number of IRs is not linear (despite
% being calculated like that) since some parts of the algorithm will be
% parallelized by Matlab already.

end
