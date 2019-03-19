% ons = AKonsetDetect(x, US, onset_threshold, thresh_mode, applyLowpass)
% Simple threshold based onset detection on up-sampled impulse responses.
% Can for example be usefull when estimating the broad band ITD.
% See AKplotDemo.m and AKp.m for examples
%
% INPUT:
% x               - data [samples x channels]
% US              - upsampling faktor US = {2,3,4...} (Default = 10).
%                   Passing 0, 1, or false will omit upsampling.
% onset_threshold - threshold for onset detection in dB (Default = -6)
% thresh_mode     - 'rel': for each channel of x, the threshold is defined
%                          by the first sample that is less than
%                          -abs(onset_threshold) dB below the max of that
%                          channel
%                   'abs': first the max values of all channels are found.
%                          Second the smalles of these values is taken,
%                          termned 'min_max'.
%                          Third the onset is defined as the first sample
%                          of each channel of x that is higher than:
%                          min_max - abs(onset_threshold)
%                   (Default = 'rel')
% applyLowpass    - two element vector [f_c fs] that gives frequency in Hz
%                   at which an 8th order Butterworth lowpass is applied to
%                   x before onset detection, and the sampling frequency.
%                   Pass false to disable (default)
%
% OUTPUT
% ons             - detected onsets in samples referring to orriginal
%                   sample rate before upsampling
%
% fabian.brinkmann@tu-berlin.de
% Audio Communication Group, TU Berlin
%
% NOTE:
% using thresh_mode = 'rel' might lead to errors there are reflections that
% exceed the level of the direct sound. In this case either 
% use thresh_mode = 'abs' or pass only the part of the impulse responses
% that contains the direct sound.
%
% TBD:
% A third thresh_mode could try to autmatically detect the level of the
% noise floor and set the onset_threshold at some point between the noise
% floor and min_max
% Onset detection using the cross-correlation between min-phase and
% original IR could be implemented. This is robust for anechoic IRs, but
% might fail elsewise

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
function ons = AKonsetDetect(x, US, onset_threshold, threshold_mode, applyLowpass)

if ~exist('US', 'var')
    US = 10;
elseif ~US
    US = 1;
end
if ~exist('onset_threshold', 'var')
    onset_threshold = 6;
end
if ~exist('threshold_mode', 'var')
    threshold_mode = 'rel';
end
if ~exist('applyLowpass', 'var')
    applyLowpass = false;
end

if strcmpi(threshold_mode, 'abs')
    onset_threshold = db(min(max(abs(x)))) - abs(onset_threshold);
end

% low-pass
if any(applyLowpass)
    % design lowpass
    Hd = fdesign.lowpass('N,F3dB', 8, applyLowpass(1)/applyLowpass(2)*2);
    Hd = design(Hd, 'butter');
    
    % check stability
    if ~isstable(Hd)
        error('AKonsetDetec:Input', 'Lowpass filter is not stable. Consider a difference cut off frequency.')
    end
    
    % apply lowpass
    x = filter(Hd, x);
end

ons = zeros(1, size(x, 2));

for r=1:size(x,2)
    % upsampling
    if US > 1
        x_us=interp(x(:,r),US);
    else
        x_us = x(:,r);
    end
    
    if strcmpi(threshold_mode, 'rel')
        % Maximum of x
        [MAX, pos] = max(abs(x_us));
        % find position were where onset threshold is reached
        ons(r) = find(abs(x_us(1:pos)) >= MAX*10^(-abs(onset_threshold)/20), 1, 'first');
    elseif strcmpi(threshold_mode, 'abs')
        % find position were where onset threshold is reached
        ons(r) = find(abs(x_us) >= 10^(onset_threshold/20), 1, 'first');
    else
        error([threshold_mode ' is not a valid for input argument ''threshold_mode'''])
    end
end

ons = ons/US + 1;
