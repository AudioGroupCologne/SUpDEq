% sm = AKfractOctSmooth(x,[operation,fs,frac,filterlen,mode,band,b,do_plot]);
%
% Applies fractional octave band smoothing on multichannel spectrum 
% (delivered as single-sided spectrums, complex or absolute values 
% are possible, complex ones will be converted to absolute values)
%
% To smooth a noisy spectrum try this
% x = fft(AKnoise)
% X = AKboth2singleSidedSpectrum(x)
% Y = AKfractOctSmooth(X,'welti',fs,frac)
%
%
%  INOUT operation - 'welti'(default) 
%
%     output:
%        sm        - smoothed single-sided spectrum with the
%                    same length as the input spectrum
%     input:
%        x         - single or multi channel single sided spectra to be
%                    smoothed (complex or absolute values allowed,
%                    1 channel = 1 column)      
%        fs        - sampling frequency, default = 44100
%        frac      - fraction of octave for smoothing, default = 3 (third octave smoothing)
%        do_plot   - [0/1], default = 0, plots only first channel
%
%                    ----------------------
%
%  INPUT operation - 'schaerer'    
%                  - takes multichannel impulse responses
%                    (timesignals) as input, windows them with a filter
%                    of length 'filterlen' and outputs a single-sided complex
%                    smoothed spectrum with half the length of the
%                    filter window
%                  - very slow, especially for big window lengths (>2048samples)
%
%  O U T P U T:
%        sm        - smoothed complex single-sided spectrum with half the length of the
%                    window declared in 'filterlen'
%     input:
%        x         - multichannel impulse response (orientation?)
%        fs        - samplerate
%        frac      - fraction of octave in case octave-width smoothing ischoosen
%        filterlen - length of the window applied on the impulse response
%        mode      - 'cpm' for complex smoothing (phase and amplitude) or 'amp'
%                    for smoothing amplitude only
%        band      - 'oct' for octave, 'bark' for bark or 'ERB'
%                    for Equivalent Rectangular Bandwidth
%        b         - coefficient for smoothing function
%                    W = b - (1-b)*cos(2*pi*k/N), for k = [0,filterlen-1]. The
%                    default value of 0.5 produces a Hann window function
%
%-------------------------------------------------------------------------%
% Additional information to operational mode 'schaerer':
%
% Calculates the smoothed transfer function to an input [h]. Returns the
% smoothed transfer function H_smooth.
% The smoothed function is calculated by convolving the transfer function
% of h with a window funktion W. The width of the window increases with
% increasing frequency, [band] determines the frequency scale used (fraction
% of octave, bark or ERB). The convolution is either conducted with the
% magnitude transfer function of h or the complex transfer function. In the
% latter case, the window is also applied to the time domain, increasing in
% width with decreasing frequency. In the case of equivalent complex
% smoothing, H_smooth is a combination of the amplitude smoothed magnitude
% and the complex smoothed phase.
% [Hatziantoniou PD, Mourjopoulos JN (2000) Generalized Fractional-Octave
% Smoothing of Audio and Acoustic Responses]
%-------------------------------------------------------------------------%
%
% ISSUES: - misses defaults for mode 'schaerer'
%
% (C)Andreas Rotter, geeprombolo@gmail.com,
%    TU Berlin, audio communication group, 2009
%    using functions from Todd Welti and Andre Giese (see below)
%-------------------------------------------------------------------------%

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

function sm = AKfractOctSmooth(x,operation,fs,frac,filterlen,mode,band,b,do_plot)


if nargin == 1
    operation = 'welti';
end

if strcmp(operation,'welti') == 1

    % check number of function args, set defaults
    narginchk(2,9)

    if nargin == 2
        fs = 44100;
        frac=3;
        do_plot=0;
    end

    if nargin == 3
        frac=3;
        do_plot=0;
    end

    if nargin == 4 || nargin == 5 || nargin == 6 || nargin == 7 || nargin == 8
        do_plot=0;
    end

    % amount of channels included in input m
    [len, wid]=size(x);

    % make input row-wise absolute, if its complex
    if ~isreal(x)
        x=abs(x);
    end

    
    % skip smoothing
    if frac == 0
        sm = x;
        return;
    end
    
    % initialize further parameters for smoothing:

    % frequency support vector, depends on length of the input vector and sample rate
    f=linspace(0,fs/2,len);

    % oct length > 1 is required, length depends on the amount of transfreqs
    % (usually only one transfreq required)
    oct = 1/frac;

    % smoothing, for each channel
    sm = zeros(size(x));
    for i=1:wid
        sm(:,i)=smoothnew3(x(:,i),oct,f(end),f);
    end

    % plotting
    if do_plot
        figure;
        semilogx(f,20*log10(x(:,1)),'k');
        hold on
        semilogx(f,20*log10(sm(:,1)),'r','LineWidth',2)
        grid on
        title('... only channel 1 ...');
        xlabel('frequency [Hz]')
        ylabel('magnitude [dB]')
        % background white
        set(gcf,'Color',[1 1 1])
    end

    %-------------------------------------------------------------------------%

elseif strcmp(operation,'schaerer') == 1
    
    if size(x,2)>size(x,1)
        x=x';
    end
    
    setup.fs=fs;
    setup.channels=size(x,2);
    
    narginchk(7, 9)
    
    if  (strcmp(band,'ERB') == 1) || (strcmp(band,'bark') == 1)
        
        b = frac;
        
    end
    
    if strcmp(band,'oct') == 1
        
        if nargin == 7
            b = 0.5;
        end
    end
    
    %-------------------------------------------------------------------------
    % N = length(x); <- this code was unused and commented, FB
    
    % n-point fft of h, with n=filterlen
    % h = pre_windowing(x,filterlen); <- this code was unused and commented, FB
    H = fft(x,filterlen);
    frac=1/frac;
    
    if strcmp(mode,'amp') == 1
        H2 = abs(H).^2;
        disp('Amplitude smoothing')
    end
    if strcmp(mode,'cmp') == 1
        H_r = real(H);
        H_im = imag(H);
        disp('Complex smoothing')
    end
    if strcmp(mode,'eqcmp') == 1
        H_r = real(H);
        H_im = imag(H);
        H2 = abs(H).^2;
        disp('Equivalent complex smoothing')
    end
    
    %----------------------------------
    %smoothing
    
    % generate vector with bandwidths according to the selected scale
    d_f = setup.fs/filterlen;
    
    switch band
        
        case 'oct'
            disp(strcat('bandwidth: ',num2str(frac),' octave'))
            bandwidth = zeros(filterlen/2+1,1);
            for k = 1:filterlen/2+1
                f_up = 10^(3/10)^(0.5*frac);
                f_low = 10^(-3/10)^(0.5*frac);
                bandwidth(k)= k*d_f*(f_up - f_low);
            end
            
        case 'bark'
            disp('bandwidth: critical band')
            bandwidth = zeros(filterlen/2+1,1);
            for k = 1:filterlen/2+1
                bandwidth(k) = 25 + 75*(1 + 1.4*(k*d_f/1000)^2)^0.69;
            end
            
        case 'ERB'
            disp('bandwidth: ERB')
            bandwidth = zeros(filterlen/2+1,1);
            for k = 1:filterlen/2+1
                bandwidth(k) = 24.7*(4.37*(k*d_f/1000) + 1);
            end
    end
    
    m = zeros(filterlen/2+1,1);
    for k = 1:filterlen/2+1
        m(k) = 0.5*bandwidth(k)/d_f;
    end
    m = ceil(m);
    clear bandwidth
    
    % generate windowfunction W
    W = zeros(length(m),filterlen);
    for j = 1:length(m)
        for k = 0:m(j)
            W(j,k+1) = (b - (1-b)*cos((pi/m(j))*(m(j)-k))) ./ (2*b*(m(j)+1) - 1);
        end
        for k = (filterlen-m(j)):(filterlen-1)
            W(j,k+1) = (b - (1-b)*cos((pi/m(j))*((filterlen-m(j))-k))) ./ (2*b*(m(j)+1) - 1);
        end
    end
    
    % generate circular convolution matrix H_circ
    if strcmp(mode,'cmp') == 1 || strcmp(mode,'eqcmp') == 1
        H_circ_r = zeros(filterlen,filterlen,setup.channels);
        H_circ_im = H_circ_r;
        for j = 1:setup.channels
            for l = 0:filterlen-1
                H_circ_r(l+1,:,j) = circshift(H_r(:,j),l);
                H_circ_im(l+1,:,j) = circshift(H_im(:,j),l);
            end
        end
    end
    if strcmp(mode,'amp') == 1 || strcmp(mode,'eqcmp') == 1
        H_circ2 = zeros(filterlen,filterlen,setup.channels);
        for j = 1:setup.channels
            for l = 0:filterlen-1
                H_circ2(l+1,:,j) = circshift(H2(:,j),l);
            end
        end
    end
    clear H_r H_im H_r2 H_im2
    
    % circular convolution of W and H_circ
    if strcmp(mode,'cmp') == 1 || strcmp(mode,'eqcmp') == 1
        H_smooth_circ_r = zeros(filterlen/2+1,filterlen,setup.channels);
        H_smooth_circ_im = H_smooth_circ_r;
        for j =1:setup.channels
            H_smooth_circ_r(:,:,j) = W * H_circ_r(:,:,j);
            H_smooth_circ_im(:,:,j) = W * H_circ_im(:,:,j);
        end
    end
    if strcmp(mode,'amp') == 1 || strcmp(mode,'eqcmp') == 1
        H_smooth_circ2 = zeros(filterlen/2+1,filterlen,setup.channels);
        for j =1:setup.channels
            H_smooth_circ2(:,:,j) = W * H_circ2(:,:,j);
        end
    end
    clear H_circ_r H_circ_im H_circ_r2 H_circ_im2 W
    
    % derive the smoothed function from H_smooth_circ
    U = zeros(filterlen,filterlen);
    v = zeros(1,filterlen/2+1);
    if strcmp(mode,'cmp') == 1 || strcmp(mode,'eqcmp') == 1
        H_smooth_r = zeros(filterlen,setup.channels);
        H_smooth_im = H_smooth_r;
        for l = 1:setup.channels
            for k = 1:filterlen/2+1
                U(k,k) = 1;
                v(m(k)) = 1;
                H_smooth_r(k,l) = sum(v*H_smooth_circ_r(:,:,l)*U);
                H_smooth_im(k,l) = sum(v*H_smooth_circ_im(:,:,l)*U);
                U(k,k) = 0;
                v(m(k)) = 0;
            end
        end
    end
    if strcmp(mode,'amp') == 1 || strcmp(mode,'eqcmp') == 1
        H_smooth2 = zeros(filterlen,setup.channels);
        for l = 1:setup.channels
            for k = 1:filterlen/2+1
                U(k,k) = 1;
                v(m(k)) = 1;
                H_smooth2(k,l) = sum(v*H_smooth_circ2(:,:,l)*U);
                U(k,k) = 0;
                v(m(k)) = 0;
            end
        end
    end
    clear U v H_smooth_circ_r H_smooth_circ_im
    
    if strcmp(mode,'cmp') == 1 || strcmp(mode,'eqcmp') == 1
        for l = 1:setup.channels
            for k = filterlen/2+2:filterlen
                H_smooth_r(k,l) = H_smooth_r(filterlen-k+2,l);
                H_smooth_im(k,l) = -H_smooth_im(filterlen-k+2,l);
            end
        end
        H_smoothc = H_smooth_r + sqrt(-1).*H_smooth_im;
    end
    if strcmp(mode,'amp') == 1 || strcmp(mode,'eqcmp') == 1
        for l = 1:setup.channels
            for k = filterlen/2+2:filterlen
                H_smooth2(k,l) = H_smooth2(filterlen-k+2,l);
            end
        end
    end
    
    switch mode
        case 'amp'
            H_smooth = sqrt(H_smooth2).*exp(sqrt(-1).*angle(H));
            
        case 'cmp'
            H_smooth = H_smoothc;
            
        case 'eqcmp'
            H_smooth = zeros(filterlen,setup.channels);
            for k = 1:setup.channels
                H_smooth(:,k) = sqrt(H_smooth2(:,k)).*exp(sqrt(-1).*angle(H_smoothc(:,k)));
            end
    end
    sm = H_smooth(1:length(H_smooth)/2,:);
else
    disp('You have to choose at least one operation mode!');
end

end

%-------------------------------------------------------------------------%
function[Hs] = smoothnew3(H,oct,stopf,f)
% ? by Andre Giese

if exist('f','var') %frequency argument vector
    df = f(2);
    %turn Hz - values to indexes
    stopf = round(stopf/df)+1;
end

% cast H to double to avoid error in smoothnew2 (F. Brinkmann 7/2012)
H = cast(H, 'double');

Hs=smoothnew2(H,oct,stopf);

%transpose if needed to get same output as input (column/row)
if size(H,1)~=stopf && size(H,2)~=1
    Hs = Hs';
end

end

%-------------------------------------------------------------------------%
function[H] = smoothnew2(H,oct,stopf)
% ? by Todd Welti

if nargin == 2
    stopf = length(H);
end

startf = 1;   % fixed at 1 for now, doesnt work right if not starting at 1
N = stopf;    % Number of samples for log warping.  This seems to work well.

spacing = 10^(log10(stopf-startf)/N);   % this is a multiplicitive factor, not an even spacing
Noct    = oct*log10(2)/log10(spacing);  % number of bins per xth octave

%log warp
logsamp = logspace(log10(startf),log10(stopf),N);
Hinterp = interp1(1:length(H),H, logsamp,'spline');

Noct_even = round(Noct/2)*2;    %so can divide by two later and use length(Noct_even) as argument to functions without uncertainty of rounding

%ADDED 11/1/08 to correct.  -3dB points of window should = Noct_even, not entire window.
% W = gausswin(Noct_even);
W = gausswin(Noct_even*2);

%extend the function to be smoothed to minimze errors at boundaries
lead = ones(1,length(W)).*Hinterp(1);
lag = ones(1,length(W)).*Hinterp(end);
Hinterp_extrap = [lead  Hinterp lag  ];

%convolve with smothing window
Htemp = fftfilt(W,Hinterp_extrap,length(W)*2)./sum(W);

Htemp = Htemp(length(lead)+.5*length(W)+1:length(lead)+.5*length(W)+length(Hinterp));   %get rid of extra bins from adding lead and lag
H = interp1(logsamp,Htemp,1:stopf,'spline');   %back to linear domain

end

% the last to functions were unused and commented out, FB

% function h_w = pre_windowing(h,len)
% 
% channels = size(h,2);
% 
% % windowing
% w = blackmanharris(len);
% w(1:len/2) = 1;
% h_w = zeros(len,channels);
% for i = 1:channels
%     h_w(:,i) = h(1:len,i).*w;
% end
% h_w = DCext(h_w,w);     % dc-extinction
% 
% end

% function h = DCext(h,w)
% 
% N = length(h);
% channels = size(h,2);
% 
% w_sum = sum(w);
% h_sum = sum(h);
% div = h_sum./w_sum;
% mult = zeros(N,channels);
% for i = 1:channels
%     mult(:,i) = div(i).*w;
%     h(:,i) = h(:,i) - mult(:,i);
% end
% 
% end
