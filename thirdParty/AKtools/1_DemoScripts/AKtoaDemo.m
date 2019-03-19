% This demo shows how to detect and remove the time of arrival (TOA) from
% impulse responses with sub-sample accurace. This might be helpfull in
% many cases, e.g.
% - for averaging or interpolating impulse responses measured at
%   different posisionts in a room
% - for removing the TOA from binaral room impulse responses, or
%   head-related impulse responses
% - for estimating the interaural time delay (ITD) in binaural impulse
%   responses
%
% The demo includes the following steps:
% 1. estimate the TOAs
% 2. fractional delaying (general)
% 3. remove TOAs
% 4. average impulse responses
%
% 12/2016 - fabian.brinkmann@tu-berlin.de

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

% in this example we take head-related impulse responses in the horizontal
% plane (left ear only)
h = AKhrirInterpolation(0:2:360, 0, 0, 'measured_ir');


%% --------------------------------------------------- 1. estimate the TOAs
%  AKtools can estimate the TOAs using onset detection or
%  cross-correlation. In the first case the first sample in the impulse
%  which exceeds a certain threshold is used, in the second case the
%  maximum cross correlation between the minimum-phase, and unaltered
%  impulse response is taken.
%  If needed, both methods can be used on the up-sampled imuplse responses
%  to achieve sub-sample accuracy. This might be helpfull when estimating
%  the ITD, because the sensiticity to humans for ITD is smaller then the
%  distance between two samples at sampling rates of 44.1 and 48 kHz

%  Set the TOA method, up-sampling factor, and smoothing
toaMethod = 'onset';   % 'onset', or 'cc'. The first method is faster the 
                       % second method might be more robust in some cases
upSample  = 10;        % 10 times up-sampling results in a quantization of
                       % 2.3 micro seconds at 44.1 kHz sampling rate.
smoothing = 1;         % sometimes outlier occur in the TOA estimation.
                       % These can be smoothed out if necessary, whereby a
                       % value of 1 denotes no smoothing and a value of 0
                       % heave smoothing - smoothing is not needed for this
                       % example.
doPlot    = true;      % the plot will show the estimated TOAs plotted on
                       % top of the impulse responses passed to AKtoa

% estimate the TOAs in samples with default parameters - for adjusting the
% parameters see the help of AKtoa.m
toa = AKtoa(h, toaMethod, false, upSample, smoothing, doPlot);

%% --------------------------------------- 2. fractional delaying (general)
%  In this demo we remove the TOAs using fractional delays realized using
%  Kaiser windowed sinc filters. We will have a look at the filters first
%  before we apply them to our impulse responses

% The filter quality depends on the filter order and side loba attenuation
% of the kaiser window
N = 30;     % The corresponding filters will be of length N+1
A = 60;     % side lobe attenuation in dB

% To illustrate the filters we will delay a simple pulse in steps of 0.1
% samples
x = AKdirac(100, 9, 50);
delays = .1:.1:.9;

% delay the pulse
xDelayed = AKfractionalDelay(x, delays, N, A);

% the plot shows that magnitude and group delay distortion for the default
% filter settings are negligible for frequencies roughly below 20 kHz
% (black line shows original pulse, colored lines the delayed version)
AKf(30,10)
subplot(1,2,1)
    AKp(xDelayed, 'm2d')
    AKp(x, 'm2d', 'c', 'k', 'dr', [-.5 .5], 'lw', 2)
subplot(1,2,2)
    AKp(xDelayed, 'gd2d', 'du', 'n')
    AKp(x, 'gd2d', 'c', 'k', 'dr', [49.5 51.5],'lw', 2, 'du', 'n')

%% --------------------------------------------------------- 3. remove TOAs
%  we now use the fractional delays to remove the TOAs from the
%  head-related impulse responses. For safety, we keep 10 samples before
%  the estimated onsets

hTOA = AKfractionalDelay(h, -toa+20);

% The plot shows an example of a head-related impulse response before and
% after fractional delaying
AKp(h(:,1), {'tc2d' 'ms2d'}, 'c', 'k')
AKp(hTOA(:,1), {'tc2d' 'ms2d'}, 'c', 'r')
legend('oiginal HRIR', ['delayed HRIR (' num2str(-toa(1)+20) ' samples)'], 'location', 'SouthWest')

%% ------------------------------------------- 4. average impulse responses
%  to illustrate the benefit of the TOA removal, we average the
%  head-related impulse responses.

hMean    = AKaverage(h, 'complex');
hMeanTOA = AKaverage(hTOA, 'complex');

% the plot shows massive loss in high frequencies, if averaging across the
% original impulse responses
AKp(hMean, '2a', 'c', 'k')
AKp(hMeanTOA, '2a', 'c', 'r')
legend('hMean', 'hMeanTOA', 'location', 'SouthWest')


