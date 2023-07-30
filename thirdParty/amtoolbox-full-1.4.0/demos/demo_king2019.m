%DEMO_KING2019 Demo of the auditory model of King, Varnet, and Lorenzi (2019)
%  
%   This script estimates the internal representations of a pure tone with 
%   and without a sinusoidal amplitude modulation, and plots two modulation
%   frequency bands for the band centred at 5 kHz.
%   The simulations are run using the model default parameters, several of 
%   them defined in arg_king2019.m .
%
%   Figure 1: Model output at 5 kHz
%
%   Figure 2: Bandpass-filtered output at 5 kHz
%
%   See also: king2019
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_king2019.php


%   #Author: Clara Hollomey (2021): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

display_level = 'no_debug'; % set to 'debug' to see more information, set to 'no_debug' to have less mess on your display
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
[outsig_target,fc,mfc,outs_step] = king2019(insig_target,fs,'basef',basef,'dboffset',dBFS,display_level);
outsig_ref = king2019(insig_ref,fs,'basef',basef,'dboffset',dBFS,display_level);

idx_fc = 3;
idx_mfc = find(mfc>=fmod,1,'first');
subfs = outs_step.subfs;

t = (1:size(outsig_target,1))/subfs;

YL = [-1.5 2]*1e-3;

figure;
plot(t,outsig_target(:,idx_fc,idx_mfc),'b'); hold on; grid on
plot(t,outsig_ref(:,idx_fc,idx_mfc),'r');
xlabel('Time [s]');
ylabel('Amplitude [arbitrary units]');
title(sprintf('Model output at fc=%.1f Hz, mod filter mfc=%.1f Hz\n(envelope)',fc(idx_fc),mfc(idx_mfc)))
%ylim(YL)

idx_mfc = idx_mfc-1;
figure;
plot(t,outsig_target(:,idx_fc,idx_mfc),'b'); hold on; grid on
plot(t,outsig_ref(:,idx_fc,idx_mfc),'r');
xlabel('Time [s]');
ylabel('Amplitude [arbitrary units]');
title(sprintf('Model output at fc=%.1f Hz, mod filter mfc=%.1f Hz \n (bandpass filtered)',fc(idx_fc),mfc(idx_mfc)))
%ylim(YL)


