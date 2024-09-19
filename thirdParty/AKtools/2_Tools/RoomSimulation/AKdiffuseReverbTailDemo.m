% This script demonstrates the use of AKdiffuseReverbTail.m. It generates
% colored, and decaying guassian noise that can serve for the simulation of
% late reverberation of a room.
%
% [1] James A Moorer (1979): "About this reverberation business." Computer
%     Music Journal, 3(2):13-28.
% [2] Christian Borss and Rainer Martin (2009): "An Improved Parametric
%     Model for Perception-Based Design of Virtual Acoustics." In: Proc.
%     of the 35th International AES Conference: Audio for Games. London
%
% 2014/11 fabian.brinkmann@tu-berlin.de

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

% reverberation time typical for a medium lecture room
T = [1.05 1.01 1   .95  .87  .7  .45];
% corresponding frequencies in Hz
f = [64   250  500 1000 2000 4000 8000];


%% ---------------------------------- demo using non-stationary combination
% This is the default method, that is the most accurate but also the least
% efficient. It uses a time-variant filter to apply the frequency and time-
% dependent decay.

% a) binaural corelation, frequency dependend
non_stat_a = AKdiffuseReverbTail(T, f);
% b) decorrelated, frequency dependend
non_stat_b = AKdiffuseReverbTail(T, f, 'rev_channel', 2, 'bin_coherence', false);
% c) mono, frequency dependend
non_stat_c = [non_stat_b(:,1) non_stat_b(:,1)];
% d) binaural corelation, frequency independend
non_stat_d = AKdiffuseReverbTail(T(4), f(4));

% a) binaural corelation, frequency dependend
% soundsc(non_stat_a, 44100)

% b) decorrelated, frequency dependend
% soundsc(non_stat_b, 44100)

% c) mono, frequency dependend
% soundsc(non_stat_c, 44100)

% d) binaural corelation, frequency independend
% soundsc(non_stat_d, 44100)


%% ------------------------------------ demo using fft with overlap and add
% This approach applies the time- and frequency dependent decay in over-
% lapping FFT blocks. It is reasonably accurate and fast.

% a) binaural corelation, frequency dependend
fft_ola_a = AKdiffuseReverbTail(T, f, 'decay_mode', 'fft_ola');
% % b) decorrelated, frequency dependend
fft_ola_b = AKdiffuseReverbTail(T, f, 'rev_channel', 2, 'bin_coherence', false, 'decay_mode', 'fft_ola');
% % % c) mono, frequency dependend
fft_ola_c = [fft_ola_b(:,1) fft_ola_b(:,1)];
% % % d) binaural corelation, frequency independend
fft_ola_d = AKdiffuseReverbTail(T(4), f(4), 'decay_mode', 'fft_ola');

% a) binaural corelation, frequency dependend
% soundsc(fft_a, 44100)

% b) decorrelated, frequency dependend
% soundsc(fft_b, 44100)

% c) mono, frequency dependend
% soundsc(fft_c, 44100)

% d) binaural corelation, frequency independend
% soundsc(fft_d, 44100)


%% --------------------------------------------------------- demo using fft
% This works as the above but with non-overlapping FFT blocks.

% a) binaural corelation, frequency dependend
fft_a = AKdiffuseReverbTail(T, f, 'decay_mode', 'fft');
% % b) decorrelated, frequency dependend
fft_b = AKdiffuseReverbTail(T, f, 'rev_channel', 2, 'bin_coherence', false, 'decay_mode', 'fft');
% % % c) mono, frequency dependend
fft_c = [fft_b(:,1) fft_b(:,1)];
% % % d) binaural corelation, frequency independend
fft_d = AKdiffuseReverbTail(T(4), f(4), 'decay_mode', 'fft');

% a) binaural corelation, frequency dependend
% soundsc(fft_a, 44100)

% b) decorrelated, frequency dependend
% soundsc(fft_b, 44100)

% c) mono, frequency dependend
% soundsc(fft_c, 44100)

% d) binaural corelation, frequency independend
% soundsc(fft_d, 44100)


%% ----------------------------------------------- demo using a filter bank
%  this is the common way of generating this type of reverb. A noise signal
%  is filterd in fractional octave bands, a decay curve is applied to each
%  partial signal, and all signals are summed at the end.

% a) binaural corelation, frequency dependend
filter_bank_a = AKdiffuseReverbTail(T, f, 'decay_mode', 'filt_bank');
% b) decorrelated, frequency dependend
filter_bank_b = AKdiffuseReverbTail(T, f, 'rev_channel', 2, 'bin_coherence', false, 'decay_mode', 'filt_bank');
% % c) mono, frequency dependend
filter_bank_c = [filter_bank_b(:,1) filter_bank_b(:,1)];
% % d) binaural corelation, frequency independend
filter_bank_d = AKdiffuseReverbTail(T(4), f(4), 'decay_mode', 'filt_bank');

% a) binaural corelation, frequency dependend
% soundsc(filter_bank_a, 44100)

% b) decorrelated, frequency dependend
% soundsc(filter_bank_b, 44100)

% c) mono, frequency dependend
% soundsc(filter_bank_c, 44100)

% d) binaural corelation, frequency independend
% soundsc(filter_bank_d, 44100)
