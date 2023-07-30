function outsig = sig_notchednoise(fc,fs,dur,L,bw,delta)
%SIG_NOTCHEDNOISE  Generates a notched-noise-type masker
%   Usage: outsig = sig_notchednoise(fs,fc,dur,L,bw,delta);
%          outsig = sig_notchednoise(fs,fc,dur,L,bw,[deltaL deltaR]);
% 
%   outsig = SIG_NOTCHEDNOISE(fs,fc,dur,L,bw,delta) generates a notched-noise
%   masker with duration dur (in sec) and overall level L (in dB SPL)
%   with a sampling rate of fs Hz. The deviation from center frequency
%   fc is symmetric and is given by delta such that the stopband is
%   [fc-delta*fc fc+delta*fc]. The left and right noise bands have a
%   bandwidth of bw*fc in Hz. If delta=0 then a broadband noise is
%   returned.
%
%   outsig = SIG_NOTCHEDNOISE(fs,fc,dur,L,bw,[deltaL deltaR]) generates a
%   notched-noise masker with an asymmetric configuration. deltaL and
%   deltaR denote the left and right deviations from fc, respectively.
%   In this case the stopband is [fc-deltaL*fc fc+deltaR*fc].
%
%   This notched-noise-type masker was used in psychoacoustical studies
%   investigating the auditory filters' shape (the original method is
%   described in Patterson, 1974). The noise is composed of two noise bands
%   of width bw and a stopband centered at fc with a deviation from fc
%   given by delta.
%
%   Examples:
%   ---------
%
%   The following shows the spectrum and a spectogram of a typical notched
%   noise masker used in the Patterson study:
%
%     fc    =  4000;  
%     fs    = 16000;
%     dur   =     1;
%     L     =   100;  
%     bw    =    .5;
%     delta =    .2;
%     outsig = sig_notchednoise(fc,fs,dur,L,bw,delta);
%
%     figure(1);
%     plotfftreal(fftreal(outsig),fs,100);
%
%     figure(2);
%     sgram(outsig,fs,80);
% 
%   References:
%     B. Moore, R. Peters, and B. Glasberg. Auditory filter shapes at low
%     center frequencies. J. Acoust. Soc. Am., 88(1):132--140, 1990.
%     
%     R. Patterson. Auditory filter shape. J. Acoust. Soc. Am., 55:802--809,
%     1974.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_notchednoise.php


%   #Author: Peter L. Sondergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin<6
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(delta) || ~isempty(find(delta<0,1))
  error('%s: delta must be a positive scalar or array.',upper(mfilename));
end;

%% Generate broadband Gaussian noise
% Make sure the length is even
n = round(dur*fs/2)*2;

noise = randn(n,1);

n_ramp = round(10E-3*fs);% Default = 10-ms Hanning on/off ramps

noise = rampsignal(noise,n_ramp);

%ramp = hann(n_ramp);
% Apply temporal windowing
%noise(1:n_ramp/2) = noise(1:n_ramp/2).*ramp(1:end/2)';% Onset
%noise(end-n_ramp/2+1:end) = noise(end-n_ramp/2+1:end).*ramp(end/2+1:end)';% Offset
% Zero padding to account for FIR delay (l = IR length)

l = 1024;% Length of filter impulse response
noiseZP = [zeros(l,1); noise; zeros(l,1)];

% Set overall level in dB SPL
noiseZP = scaletodbspl(noiseZP,L);    
    
%% If delta = 0 then filter is not required
if isscalar(delta) && delta == 0
    outsig = noiseZP;
end
%% Multiband filter design (FIR filter)
if isscalar(delta) && delta > 0
%     Symmetric notch
    b1l = ((fc*(1-delta-bw))*2)/fs;% Low edge of left noise band
    if b1l < 0
        b1l = 0;
    end
    b1h = ((fc*(1-delta))*2)/fs;% High edge of left noise band
    b2l = ((fc*(1+delta))*2)/fs;% Low edge of right noise band
    b2h = ((fc*(1+delta+bw))*2)/fs;% High edge of right noise band
    if b2h > 1
        b2h = 1;
    end
    % Compute and analyze filter
    f = [0 b1l b1l b1h b1h b2l b2l b2h b2h 1];
    m = [0 0 1 1 0 0 1 1 0 0];
    b = firls(l,f,m);
    % Plot for verification (uncomment if needed)
    % [h,w] = freqz(b,1,1024);
    % figure
    % plot(f,m,w/pi,abs(h),'--'),grid on
    % xlabel('Normalized frequency'), ylabel('Magnitude'), title('Filter response')
    % Apply filter:
    outsig = filter(b,1,noiseZP);
elseif ~isscalar(delta) 
    if length(delta) > 2
        error('%s: Stopband is not correctly specified.',upper(mfilename));
    end
%     Asymmetric notch
    b1l = ((fc*(1-delta(1)-bw))*2)/fs;% Low edge of left noise band
    if b1l < 0
        b1l = 0;
    end
    b1h = ((fc*(1-delta(1)))*2)/fs;% High edge of left noise band
    b2l = ((fc*(1+delta(2)))*2)/fs;% Low edge of right noise band
    b2h = ((fc*(1+delta(2)+bw))*2)/fs;% High edge of right noise band
    if b2h > 1
        b2h = 1;
    end
    % Compute and analyze filter
    f = [0 b1l b1l b1h b1h b2l b2l b2h b2h 1];
    m = [0 0 1 1 0 0 1 1 0 0];
    b = firls(l,f,m);
    % Plot for verification (uncomment if needed)
    % [h,w] = freqz(b,1,1024);
    % figure
    % plot(f,m,w/pi,abs(h),'--'),grid on
    % xlabel('Normalized frequency'), ylabel('Magnitude'), title('Filter response')
    % Apply filter:
    outsig = filter(b,1,noiseZP);
end

%% Plot results (uncomment if needed)
% fft_noise1 = fft(noiseZP)./(n+2*l);
% fft_noise2 = fft(outsig)./(n+2*l);
% figure
% % Time domain
% subplot(2,1,1)
% plot(noiseZP), hold on
% plot(outsig,'--r'), hold off
% legend('Input','Outsig')
% xlabel('Samples'),ylabel('Amplitude'),title('Time domain')
% % Frequency domain
% subplot(2,1,2)
% plot(linspace(0,fs,n+2*l),20*log10(abs(fft_noise1))), hold on
% plot(linspace(0,fs,n+2*l),20*log10(abs(fft_noise2)),'--r'), hold off
% legend('Broadband noise','Filtered noise')
% xlabel('Frequency (Hz)'),ylabel('Squared modulus (dB SPL)'),title('Frequency domain')

% eof


