function s = sig_linsweep(fs, N, range)
%SIG_LINSWEEP Create a linear frequency sweep with constant magnitude spectrum
%   Usage: s = sig_linsweep(fs, N, range)
%
%   S = sig_linsweep(FS, N) returns a sweep of length N with normalized 
%   amplitude and sampling rate FS with a perfectly constant magnitude 
%   spectrum. It covers all frequencies from 0 to FS/2 and runs from 
%   the first to the last sample.
% 
%   S = sig_linsweep(FS, N, RANGE) returns a sweep of length N with 
%   normalized amplitude and a sampling rate of FS and a perfectly 
%   constant magnitude spectrum. The signal covers all frequencies 
%   from 0 to FS/2. The sweep starts at sample RANGE(1) 
%   and ends at sample RANGE(2), i.e., RANGE must be an array
%   with two positive elements with 0 < RANGE(1) < RANGE(2) <= N.
%   A tail is encountered before and after the actual sweep range. 
% 
%   EXAMPLE 1: Create a perfect sweep and display its spectrogram:
% 
%     fs = 44100;                 % sampling rate of 44100 Hz
%     s = sig_linsweep(fs, 4*fs)  % 4 secs sweep
%     spectrogram(s, 256, 128, 256, fs, 'yaxis');
%     sound(s,fs);
% 
%   EXAMPLE 2: Create a perfect sweep in a range:
% 
%     fs = 44100;                 % sampling rate of 44100 Hz
%     s = sig_linsweep(fs, 4*fs, [fs+1,2*fs]); % 4 secs with sweep 
%                                              % from fs+1 to 2*fs
%     spectrogram(s, 256, 128, 256, fs, 'yaxis');
%     sound(s,fs);
% 
%   To create a continuous measurement stimulus, you have to stack as 
%   many length N periods together as you need. Considering the sampling
%   rate fs, you have to repeat the sweep period t*fs/N times to get a 
%   signal length of t seconds. Make sure t*fs/N is an integer.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_linsweep.php


%   #Author: Aulis Telle (2008)
%   #Author: Piotr Majdak (2017)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if exist('range','var')
    if length(range) ~= 2 ...
            || range(1) > range(2) ...
            || range(1) < 1 ...
            || range(2) > N
        error('RANGE must be an array with two elements 0 < RANGE(1) < RANGE(2) <= N.');
    end
else
    range = [1 N];
end



% initialization
isEven = (mod(N,2) == 0);
Nhalf = ceil(N/2);
ph = zeros(1, N);



% The group delay (in seconds) of every frequency bin is
% the index of the frequency bin divided by half of the 
% sampling frequency. This way the frequency corresponding to 
% frequency index k occurs at time index 2*k.
groupDelay = ((0:Nhalf-1) / N * range(2) + range(1) - 1) / (fs/2);

% Width of a frequency bin (in radians)
deltaOmega = 2 * pi * fs / N;

% The phase is calculated by integrating over the linear 
% group delay, i.e., the value of the group delay times the 
% current frequency divided by two (because we have a 
% triangular area under the groupDelay curve).
ph(1:Nhalf) = - groupDelay .* (deltaOmega * (0:length(groupDelay)-1))./2;

% Construct an odd-symetrical phase needed for a real valued signal 
% in time domain.
if isEven
    ph(Nhalf+1:end) = [0 -fliplr(ph(2:Nhalf))];
else
    ph(Nhalf+1:end) = -fliplr(ph(2:Nhalf));
end

% Calculate complex represenation with constant magnitude for IFFT 
c = 10*exp(i*ph);

% Calculate IFFT to gain sweep signal in time domain
s = real(ifft(c));

% Normalize signal to have a maximal amplitude of one
s = s / max(abs(s));


