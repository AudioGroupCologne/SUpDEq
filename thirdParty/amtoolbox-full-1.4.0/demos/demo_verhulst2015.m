%DEMO_VERHULST2015 Demo of the cochlear model calculating otoacoustic emissions
%  
%   This script computes and plot the otoacoustic emission for a 500 Hz sinusoid,
%   This is generated using the cochlear model described in verhulst2015.
%   In particular otoacoustic emissions are computed as the signal difference between the
%   sound pressure at the middle ear with model nonlinearities and irregularities enabled,
%   and the sound pressure at the middle ear in case of linear model.
%
%   Figure 1: Output of the cochlear model.
%
%   Figure 2: Simulated basilar membrane velocity for the model with linearities and nonlinearities.
%
%   License:
%   --------
%
%   This model is licensed under the UGent Academic License. Further usage details are provided 
%   in the UGent Academic License which can be found in the AMT directory "licences" and at 
%   <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
%
%   See also: verhulst2015
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_verhulst2015.php


%   #Author: Alejandro Osses (2020): primary implementation for the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0
%   #License: ugent

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

display_level = 'no_debug'; % set to 'debug' to see more information, set to 'no_debug' to have less mess on your display

%%% 1. Model parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irr_on=[1,0];
normaliseRms=[1 1]; % only used in the old way of calling verhulst2012.
subjectNo=rand(); 
fc_flag='all'; % 1000 cochlear sections
version_year = 2015;

%%% 1. Generating the input signals: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal parameters:
fs=48000; % Hz
dur=50e-3; % 50 ms
dt = 1/fs; % \Delta t in s
t=0:dt:dur-dt;
f0=500; % Hz, Carrier frequency of the signals
spl=[60,60]; % Level of the test signals

insig=zeros(length(t),2); % Memory allocation
insig(:,1)=sin(2*pi*f0*t);
insig(:,2)=insig(:,1);

% Calibration of the input signals:
for j = 1:length(spl)
    p0 = 2e-5;
    insig(:,j) = insig(:,j)/rms(insig(:,j));
    insig(:,j) = p0*10^(spl(j)/20.)*insig(:,j);
    
    % % Equivalent code using AMT functions:
    % dBFS = 94; % i.e., amplitude 1 = 1 Pa = 94 dB SPL re 2x10^{-5} Pa
    % insig(:,j) = scaletodbspl(insig(:,j),spl(j),dBFS);
end

%%% 3. Calling the model: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% New way to call the model:
OAE = []; V=[];
flag_oae = 'oae'; % to return otoacoustic emissions (disabled by default):
outsig = verhulst2015(insig,fs,fc_flag,'subject',subjectNo,'irr_on',irr_on,flag_oae,display_level, 'v'); % irregularities are on by default
model_prefix = 'verhulst2015';

OAE(:,1)  = outsig(1).oae; % Emissions with Zweig irregularities (irr_on = 1)
OAE(:,2)  = outsig(2).oae; % Emissions without Zweig irregularities (irr_on = 0)
CF = outsig(1).cf; % characteristic frequencies
idx = find(CF<1000,1,'first'); % Looking for one specific CF (arbitrary)
V(:,1) = outsig(1).v(:,idx); % Veloc. of the basilar membrane at CF(idx), irr_on = 1
V(:,2) = outsig(2).v(:,idx); % Veloc. of the basilar membrane at CF(idx), irr_on = 0

OtoacousticEmissionPa=OAE(:,1)-OAE(:,2);

%%% 4. Plotting the results: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_ms = t*1000;
figure; 
plot(t_ms,OtoacousticEmissionPa);
grid on;
xlabel('Time (ms)');
ylabel('Emission (Pa)');
title([model_prefix ' model: Otoacoustic emission for a 500 Hz sinsuoid']);

figure;
plot(t_ms,V(:,1),'b-'); hold on; grid on
plot(t_ms,V(:,2),'r--'); 
xlabel('Time (ms)')
ylabel('Basilar membrane velocity (m/s)')
title(sprintf('%s model: Basilar membrane velocity CF=%.1f Hz (bin number=%.0f)', ...
               model_prefix,CF(idx),idx));

legend('Irregularities on','Irregularities off')


