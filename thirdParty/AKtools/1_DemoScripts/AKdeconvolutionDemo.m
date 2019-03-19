% Demonstration of spectral deconvolution to obtain an impulse response
% from a recorded sweep signal.
%
% For a demonstration of sweep generation see:
% AKtestSignalsDemo.m
%
% For a demonstration of impulse response measurements see:
% AKmeasureDemo.m
%
% 10/2016 - initial dev. fabian.brinkmann@tu-berlin.de

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

%% 0. To demonstrate the use of the AK deconvolution we will first simulate
%  an impulse response measurement:

% our sampling rate in Hz
fs = 44100;

% this is the measurement signal - a simple logarithmic sweep
x = AKsweepFD('NFFT', 16, 'preset', 'log', 't_start', 100, 't_gap', 100, 'do_plot', false);

% this is is our impulse response that we want to measure - a band pass
h = AKfilter(AKdirac(2^16, 1, 100), 'bp', [150 16000], 0, fs, 4, 'butter');

% to make things more realistic we will add some noise
n = AKnoise(2^16) / 2^16;

% we will also add the frequency response of a microphone
m = AKfilter(AKdirac(2^10), 'peq', 8000, 3, fs, 1, 'hpl', 'tan');

% simulate the measurement by convolving the exitation signal with the
% impulse response. y is what we obtain after recording a sweep signal
% played back by a loudspeaker with a measurement microphone:
y = fftfilt(h+n, x);
y = fftfilt(m, y);

% plot the signals
AKf(30,40)
subplot(4,2,1)
    AKp(x, 't2d'); title 'Sweep time signal'
subplot(4,2,2)
    AKp(x, 'm2d'); title 'Sweep magnitude spectrum'
subplot(4,2,3)
    AKp(h, 't2d', 'x', [0 10]); title 'IR time signal'
subplot(4,2,4)
    AKp(h, 'm2d'); title 'IR magnitude spectrum'
subplot(4,2,5)
    AKp(n, 't2d'); title 'Noise time signal'
subplot(4,2,6)
    AKp(n, 'm2d'); title 'Noise magnitude spectrum'
subplot(4,2,7)
    AKp(m, 't2d', 'x', [-.1 1]); title 'Microphone time signal'
subplot(4,2,8)
    AKp(m, 'm2d'); title 'Microphone magnitude spectrum'



%% 1. Now we deconvolve everything with default parameters
h_1 = AKdeconv(y, x, fs, 'plot_deconv', true, 'plot_inv', true);

% The first plot shows the deconvolved impulse response and spectrum (h_1).
% The second plot shows the inverse spectrum of x. Since h is obtained by
% spectral division H = Y/X, and inverse fourier transform h=ifft(H) the
% term 1/X might boost the noise in h if it contains high energy close to 0
% Hz or fs/2. You can see this in the second plot.

%% 2. We can tackle this by restricting the dynamic of 1/X
h_2 = AKdeconv(y, x, fs, 'x_inv_dyn', 35, 'plot_deconv', true, 'plot_inv', true);

% if you compare the results from 1. and 2. you will see only subtle
% differences in the low frequency noise contained in h_estimated. This
% effect might be more drastic with real measurements which might make it
% pretty important to check the dynmaic of 1/X. Please note that if this
% alone does not work, AKdeconv can further restrict the bandwidth of 1/X
% using highPass, and lowPass paramters
% (enter 'help AKdeconv' in the command window for more info)

%% 3. Because the result h_2 still contains the microphone frequency
%  response, we will remove it using AKdeconv by applying 1/M. Again we can
%  restric the dynamic of 1/M if by using mic_inv_dyn.
h_3 = AKdeconv(y, x, fs, 'x_inv_dyn', 35, 'mic_comp', m, 'plot_deconv', true, 'plot_inv', true);

%% 4. In case you are measuring lots of impulse responses you might want to
%  truncate them right away to safe disk space
h_4 = AKdeconv(y, x, fs, 'x_inv_dyn', 35, 'mic_comp', m, 'h_trunc', 2^12, 'plot_deconv', true, 'plot_inv', true);

% Note that truncation might affect the low end of the frequency response.
% As a rule of thumb your IR should be 2 to 4 times the lengths of the
% lowest frequency you are interested in. In this case this would be about
% 40 Hz (fs/2^12 * 4)

%% 5. Plot all results in one plot

AKf
AKp([h h_1 h_2 h_3 [h_4; zeros(2^16-2^12,1)]], 'm2d', 'c', 'cyc')
legend('h', 'h_1', 'h_2', 'h_3', 'h_4', 'location', 'best')
