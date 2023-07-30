function sig = sig_whitenoiseburst(fs)
% SIG_WHITENOISEBURST burst of white noise
%
%   Generate the noise signal used for the binaural model to predict the
%   perceived direction; used by Wierstorf et al. (2013)
%
%   See also: wierstorf2013
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_whitenoiseburst.php


%   #Author: Hagen Wierstorf (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


noise_length = 700;
pause_length = 300;
slope_length = 20;
% calculate samples from ms
noise_samples = round(noise_length/1000*fs);
pause_samples = round(pause_length/1000*fs);
slope_samples = round(slope_length/1000*fs);
% on- and offset slope
win = hann_window(slope_samples,slope_samples,noise_samples);
% bandpass filter
n=4;      % 2nd order butterworth filter
fnq=fs/2;  % Nyquist frequency
wn=[125/fnq 20000/fnq];    % butterworth bandpass non-dimensional frequency
[b,a]=butter(n,wn); % construct the filter
%
pause_sig = zeros(round(pause_samples/2),1);
sig_noise = noise(noise_samples,1,'white') .* win;
sig = [pause_sig; sig_noise; pause_sig];
sig = filtfilt(b,a,sig);
sig = sig(1:end-round(pause_samples/2));


