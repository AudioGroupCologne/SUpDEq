% Adaptive low frequency extention (ALFE) for measured impulse responses.
% fd stands for frequency domain
%
% Low frequency extension/extrapolation is done by extrapolating magnitude
% and phase response. For details, see:
% Xie, Bosun: 'On the low frequency characteristics of head-related
% transfer function'. In: Chinese J. Acoust. 28(2):1-13 (2009).
%
%
% INPUT (default parameters)
% h                  - impulse respone single/multi channel
%                      [one channel = one column]
% fs (44100)         - sampling frequency in Hz
% NFFT (size(h,1))   - if NFFT is larger/smaller than size(h,1), h is zero
%                      padded or truncated before FFT. Output length of
%                      HRIRs will be NFFT.
% f_link ([200 400]) - frequency range in Hz used for calculating the
%                      extrapolation. HRTFs should have a constant
%                      magnitude response in this region
% blend (1)          - Apllies a fade/blend between measured and
%                      extrapolated magnitude response values
% L_target (NaN)     - Target level for extrapolation of magnitude
%                      response in dB. If NaN, the magnitude response
%                      below f_link(1) is set to the average in the range
%                      of f_link. Otherwise a linear extrapolation from
%                      the average in the range of f_link to L_target in dB
%                      is performed.
% ALFE_plot (0)      - plots resulting  signals
%
% Fabian Brinkmann, Audio Communication Group, TU Berlin, 04/2013

function h_LFE =  ALFE_fd(h, varargin)
%% ------------------------------------------ 1a. define default parameters
inputnames = {...
'fs', 'NFFT'...
'f_link', 'blend' 'L_target' ...
'ALFE_plot' ...
};

% Define default values
def.fs        = 44100;
def.NFFT      = size(h,1);

def.f_link    = [200 400];
def.blend     = 1;
def.L_target  = NaN;

def.ALFE_plot = 0;

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


% ---------------------------------------------------- 2 helpfull variables
% save unprocessed input for plot
h_org = h;

%Get single sided spectrum
H     = fft(h, NFFT);
H     = AKboth2singleSidedSpectrum(H);
H_abs = abs(H);
H_ang = unwrap(angle((H)));

%Calculate indeces of linking frequencies
n_link = freq2n(f_link, NFFT, fs);

% get frequency vector
f = n2freq((1:size(H,1))', NFFT, fs);

% ------------------------------------------- 3 calculate ALFE level from h
ALFE_abs = mean(H_abs(n_link(1):n_link(2),:));
if isnan(L_target)
    ALFE_abs = repmat(ALFE_abs, [n_link(2), 1]);
else
    tmp = ALFE_abs;
    ALFE_abs = zeros(n_link(2), numel(tmp));
    for nn = 1:numel(tmp);
        ALFE_abs(:,nn) = linspace(10^(L_target/20), tmp(nn), n_link(2))';
    end
    clear tmp
end

% ----------------------------------------------------- 4 extrapolate phase
% obtain frequncy vector
ALFE_ang = interp1(f(n_link(1):n_link(2)) , H_ang(n_link(1):n_link(2),:), f(1:n_link(2)), 'linear', 'extrap');

% --------------------------- 5 couple extrapolated values to data and ifft
H_abs(1:n_link(1)-1,:) = ALFE_abs(1:n_link(1)-1,:);
if blend
    blend = linspace(1,0,n_link(2)-n_link(1)+3)';
    blend = blend(2:end-1);
    blend = repmat(blend, [1 size(H_abs,2)]);
    
    H_abs(n_link(1):n_link(2),:) = ...
        blend     .* ALFE_abs(n_link(1):n_link(2),:) + ...
        (1-blend) .* H_abs(n_link(1):n_link(2),:);
end

H_ang(1:n_link(2),:) = ALFE_ang;

% get complex spectrum
H_LFE = H_abs .* exp(1j*H_ang);
% get both sided spectrum and impulse responses
H_LFE = AKsingle2bothSidedSpectrum(H_LFE, mod(NFFT+1,2));
h_LFE = ifft(H_LFE, 'symmetric');



% ----------------------------------------------------------------- 5 plots
if ALFE_plot 
    fID = figure;
    fWidth = 30; fHeight = 20;
    set(fID,'PaperUnits', 'centimeters');
    set(fID,'Units', 'centimeters');
    % paper size for printing
    set(fID, 'PaperSize', [fWidth fHeight]);
    % location on printed paper
    set(fID,'PaperPosition', [.1 .1 fWidth-.1 fHeight-.1]);
    % location and size on screen
    set(fID,'Position', [0 0 fWidth fHeight]);
    % set color
    set(fID, 'color', [1 1 1])
    
    % restrict limits of y-axis
    y_max = max(20*log10(abs(fft(h_org))));
    y_max = max(y_max + 5-mod(y_max,5));
    y_min = 20*log10(ALFE_abs(1));
    y_min = y_min - 40 - mod(y_min, 5);
    
    % spec fullrange
    subplot(2,2,1)
        AKp(h_org, 's2d', 'c', [.7 .7 .7])
        AKp(h_LFE, 's2d', 'c', 'k', 'x', [20 fs/2])
        ylim([y_min y_max])
        title({'Magnitude response' '(grey: input; black: output)'})
    % phase fullrange
    subplot(2,2,3)
        AKp(h_org, 'ph2d', 'c', [.7 .7 .7])
        AKp(h_LFE, 'ph2d', 'c', 'k', 'x', [20 fs/2])
        
    % restrict limits of y-axis
    y_max = 20*log10(ALFE_abs(1));
    if mod(y_max,10) > 5
        y_max = y_max + 15-mod(y_max,10);
    else
        y_max = y_max + 10-mod(y_max,10);
    end
    y_min = 20*log10(ALFE_abs(1));
    if mod(y_min, 10) < 5
        y_min = y_min - 5 - mod(y_min, 10);
    else
        y_min = y_min - mod(y_min, 10);
    end
        
    % spec crossover
    subplot(2,2,2)
        AKp(h_org, 's2d', 'c', [.7 .7 .7])
        AKp(h_LFE, 's2d', 'c', 'k', 'x', [f_link * 2^-(2/2) f_link * 2^(2/2)])
        % plot limits for averaging level and group delay and plot f_link
        hold on
        plot(n2freq([n_link(1) n_link(1)], NFFT, fs), [y_max y_min], '--g')
        plot(n2freq([n_link(2) n_link(2)], NFFT, fs), [y_max y_min], '--g')
        ylim([y_min y_max])
        title({'Magnitude response' '(green: linking frequencies)'})
    % phase crossover
    subplot(2,2,4)
        AKp(h_org, 'ph2d', 'c', [.7 .7 .7])
        AKp(h_LFE, 'ph2d', 'c', 'k', 'x', [f_link * 2^-(2/2) f_link * 2^(2/2)])
        % plot limits for averaging level and group delay and plot f_link
        hold on
        plot(n2freq([n_link(1) n_link(1)], NFFT, fs), [-200 200], '--g')
        plot(n2freq([n_link(2) n_link(2)], NFFT, fs), [-200 200], '--g')
end
end