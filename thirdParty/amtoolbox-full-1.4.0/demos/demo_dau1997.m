%DEMO_DAU1997 Demo of the auditory model of Dau, Kollmeier, and Kohrausch (1997)
%  
%   This script estimates the internal representations of a pure tone with 
%   and without a sinusoidal amplitude modulation, and plots two modulation
%   frequency bands for the band centred at 5 kHz.
%   The simulations are run using the model default parameters.
%
%
%   Figure 1: The above figure shows the model output at 5 kHz for the modulation filter centred at 10 Hz.
%
%      This time signal corresponds to the 'venelope' post-processing, 
%      which keeps the Hilbert envelope of the modulation filter (similar
%      for other modulation bands with mfc>=10 Hz).
%
%
%   Figure 2: The above figure shows the bandpass-filtered output at 5 kHz for the modulation filter centred at 5 Hz.
%      
%      This time signal corresponds to the band-passed modulation
%      signal. The same processing is adopted by the model for modulation
%      bands such that mfc>2.5 Hz and mfc<10 Hz.
%
%   See also: dau1997 demo_king2019
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_dau1997.php


%   #Author: Alejandro Osses (2023)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

fs = 48000; % Hz, sampling frequency
lvl = 70;
dBFS = 100; % 0 dBFS is equal to 100 dB

fc_carrier = 5000; % Hz, frequency of the carrier
fmod = 10; % Hz, frequency of the modulator
mdepth = 1; % fully modulated signal

dur      = 500e-3; % stimulus duration in seconds
dur_ramp = 2.5e-3; % s, duration of the ramp

N_samples = round(dur*fs);
dur_ramp_samples = round((dur_ramp)*fs);

% Creating a cosine ramp:
ramp = ones(N_samples,1);
ramp(1:dur_ramp_samples)         = rampup(dur_ramp_samples);
ramp(end-dur_ramp_samples+1:end) = rampdown(dur_ramp_samples);

% AM stimulus and calibration:
t = (0:N_samples-1)/fs; t=t(:);
carrier = sin(2*pi*fc_carrier*t); % starts at phase = 0
env = (1 + mdepth * sin(2*pi*fmod*t-pi/2) ); % modulator starts at minimum (phase=-pi/2)

insig_target = env .* carrier; % Amplitude-modulated signal
insig_target = scaletodbspl(insig_target,lvl,dBFS);
insig_target = ramp.*insig_target; % ramp applied after the calibration

insig_ref = carrier;
insig_ref = scaletodbspl(insig_ref,lvl,dBFS);
insig_ref = ramp.*insig_ref; % ramp applied after the calibration

% Model parameters:
basef = fc_carrier; % Frequency of the central auditory filter.
[outsig_target,fc,mfc] = dau1997(insig_target,fs,'basef',basef,'dboffset',dBFS);
outsig_ref = dau1997(insig_ref,fs,'basef',basef,'dboffset',dBFS);

idx_fc = find(basef == round(fc)); 
idx_mfc = find(mfc>=fmod,1,'first');

t = (1:size(outsig_target{1},1))/fs;

figure;
plot(t,outsig_target{idx_fc}(:,idx_mfc),'b'); hold on; grid on
plot(t,outsig_ref{idx_fc}(:,idx_mfc),'r');
xlabel('Time [s]');
ylabel('Amplitude [Model Units]');
title(sprintf('Model output at fc=%.1f Hz, mod filter mfc=%.1f Hz\n(envelope)',fc(idx_fc),mfc(idx_mfc)))

idx_mfc = idx_mfc-1;
figure;
plot(t,outsig_target{idx_fc}(:,idx_mfc),'b'); hold on; grid on
plot(t,outsig_ref{idx_fc}(:,idx_mfc),'r');
xlabel('Time [s]');
ylabel('Amplitude [Model Units]');
title(sprintf('Model output at fc=%.1f Hz, mod filter mfc=%.1f Hz \n (bandpass filtered)',fc(idx_fc),mfc(idx_mfc)))


