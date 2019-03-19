% Demonstration of test signal generation
%
% 1. Dirac impulses
% 2. Noise and pulsed noise
% 3. Sine signals
% 4. Sine sweeps
%
% 10/2016 - initial dev. fabian.brinkmann@tu-berlin.de
%
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

% suppress unused variables warning for demo script
%#ok<*NASGU>

% set the sampling rate in Hz
fs = 44100;


%% --------------------------------------------------------------- 1. Dirac
%  A dirac is a simple pulse. This would be two lines of Matlab code, now
%  its only one:

x = AKdirac(1024, 1, 256);                      % generate
AKp(x, 't2d', 'dr', [-.2 1.2], 'xu', 'n')      % plot


%% ---------------------------------------------- 2. Noise and pulsed noise
%  Gaussian noise with pink or white spectrum

% a) 1 second gaussian noise with pink spectrum
x = AKnoise;

% plot and listen
AKf
subplot(2,1,1)
    AKp(x, 't2d')
subplot(2,1,2)
    AKp(x, 'm2d')

% soundsc(x, fs)

%% b) 1 second gaussian noise with white spectrum
x = AKnoise(fs, 1, 'white');

% plot and listen
AKf
subplot(2,1,1)
    AKp(x, 't2d')
subplot(2,1,2)
    AKp(x, 'm2d')
    
% soundsc(x, fs)

%% c) pulsed noise:
%  - 1 second duration
%  - sampling rate 44100 Hz
%  - pulse duration 0.2 seconds
%  - pulse fade in/out 0.02 seconds
%  - pause between pulses 0.1 seconds
% final length is adjusted to not cut a sequence of noise
x = AKpulsedNoise(1, fs, .2, .02, .1);

% plot and listen
AKp(x, 't2d')
% soundsc(x, fs)


%% ---------------------------------------------------------------- 3. Sine
%  Pure tones can be generated this way

% harmonic pure tones with a fundamental of 200 Hz with a length of 1
% second and a sampling rate of 44100 Hz.
% If the last argument is 'true' the sines end with an amplitude close to
% zero and the length of x is adjusted accordingly. If the last argument is
% 'false' the sines are cut regardless of the ampitude.
x = AKsine(.5, 200:200:1000, fs, false);

% apply amplitude decay
x = AKm(x, logspace(0, -3, size(x,2)), '*');

% fade in and out
x = AKfade(x, size(x,1), .02*fs, .4*fs);

% plot
AKp(x, 't2d', 'c', 'cyc')

% listen
% soundsc(sum(x, 2), fs)


%% --------------------------------------------------------- 4. Sine sweeps
%  These are sweep signals that change their frequency across time. They are
%  commonly used to measure acoustic systems such as room impulse responses
%  or head-related impulse responses. This is due to their good properties
%  - high energy in relation to the maximum amplitude
%  - possiblity of measuring non-harmonic products
%  - fast convergence in adaptive system identification
%  - arbitrary coloration for constant snr (signal to noise ration) across 
%    in frequency
%
%  lets have a close look

%  a) full range linear sweep, designed in the time domain (TD)
%     this is the simplest sweep, however not often used
x = AKsweepTD('lin', 2^16, [0 fs/2]);

% plot and listen
AKf
subplot(2,1,1)
    AKp(x, 't2d')
subplot(2,1,2)
    AKp(x, 'm2d')
    
% soundsc(x, fs)

% NOTE: If you measure an acoustic system, you will have to record some
% time after the sweep stopped to capture the decay of the system. The time
% domain sweep generation does not account for this by itself. You have to
% do this manually.


%% b) full range exponential sweep, designed in the time domain (TD)
%     the frequency changes exponentially accross time. Also called log
%     sweep, because the group delay changes logarithmically across time.
%     This can be used for measureing non-harmonic impulse responses.
%     Moreover its coloration (increased energy at low frequencies) matches
%     the noise in many acoustic environments that usually have more noise
%     at lower frequencies. The log. sweep could thus be seen as a first
%     approximation of constant snr measurement
x = AKsweepTD('exp', 2^16, [30 fs/2]);

% plot and listen
AKf
subplot(2,1,1)
    AKp(x, 't2d')
subplot(2,1,2)
    AKp(x, 'm2d')
    
% soundsc(x, fs)

% NOTE: If you measure an acoustic system, you will have to record some
% time after the sweep stopped to capture the decay of the system. The time
% domain sweep generation does not account for this by itself. You have to
% do this manually.


%% c) full range linear sweep, designed in the frequency domain (FD)
%     This sweep is often used for adaptive measurements, because of its
%     fast convergence speed. Due to this property, it is also called
%     'perfect sweep'.
x = AKsweepFD('NFFT', 2^12, 'preset', 'perfect');

% listen:
% the sweep is cyclic and can be repeated, which is done in adaptive
% measurements

%soundsc(repmat(x, 20, 1), fs)


%% d) FD sweeps with arbitrary magnitude response for constant snr measurements
%     To obtain a constant snr across frequency, the magnitude spectrum of
%     the sweep should be adjusted according the the background noise, the
%     speaker frequency response, and the microphone frequency response.
%     This is illustrated here

% get the noise floor
% (we simply take pink noise here, in real environments you could record this)
n = AKnoise;

% get the speaker frequency response
% (we model this by a simple parametric eq, and a high-pass)
s = AKfilter(AKdirac, 'peq', 1000, -4, fs, .5, 'hpl', 'cos');
s = AKfilter(s, 'hp', 150, 0, fs, 2, 'butter');

% design the sweep
% Note that the speaker's frequency response is inverted, and that you can
% restric the dynamic range of this process with 'speaker_dynamic'. We also
% apply a high-pass in the sweep design, because our speaker starts at 150
% Hz.
x = AKsweepFD('noise_ir', n, 'speaker_ir', s, 'speaker_dynamic', 5, 'highPass', [50 4 1], 'lowshelf_apply', 0);

% listen
% soundsc(x, fs)

% NOTE: If you measure an acoustic system, you will have to record some
% time after the sweep stopped to capture the decay of the system. The
% freequency domain sweep generation accounts for this by itself. You can
% specify the silence after the sweep with 't_gap'
