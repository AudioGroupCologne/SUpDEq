%DEMO_VERHULST2018 Demo of the auditory brainstem model described by Verhulst et al (2018)
%  
%   This script computes the responses to one click of single polarity.
%
%   The following figures are obtained:
%
%   - Characteristic frequencies of equidistant points along the cochlea.
%
%   - Ear-canal pressure signals that can be used to simulate 
%     otoacoustic emissions (see also demo_verhulst2012.m).
%
%   - Simulated basilar membrane velocity (v) or inner-hair-cell
%     voltage (vihc).
%
%   - Simulated unit responses to HSR, MSR, LSR neurons, auditory
%     nerve responses, and cochlear nucleus (CN) and inferior colliculus (IC) neurons.
%
%   - Population responses waave I, III, V, and envelope following responses.
%
%
%   The current demo is based on the scripts ExampleSimulation.m and 
%   ExampleAnalysis_function.m from https://github.com/HearingTechnology/Verhulstetal2018Model
%
%   License:
%   --------
%
%   This model is licensed under the UGent Academic License. Further usage details are provided 
%   in the UGent Academic License which can be found in the AMT directory "licences" and at 
%   <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.%
%
%   Figure 1: Characteristic frequencies of equidistant points along the cochlea.
%
%   Figure 2: Ear-canal pressure signals that can be used to simulate otoacoustic emissions.
%
%   Figure 3: Simulated basilar membrane velocity (v) or inner-hair-cell voltage (vihc).
%
%   Figure 4: Simulated unit responses to HSR, MSR, LSR neurons, auditory nerve responses, and cochlear nucleus (CN) and inferior colliculus (IC) neurons.
%
%   Figure 5: Population responses waave I, III, V, and envelope following responses.
%
%   See also: verhulst2018, verhulst2012
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_verhulst2018.php


%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alejandro Osses (2020): primary implementation based on https://github.com/HearingTechnology/Verhulstetal2018Model
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%%% Generic model configuration:
fc_flag='all'; % 1000 cochlear sections
% fc_flag = 'half'; % for 500 cochlear sections
% fc_flag = 'abr';  % for 401 cochlear sections
numH = 13; % Nr. of high-spontaneous   rate neurons to be simulated (default=13)
numM =  3; % Nr. of middle-spontaneous rate (default=3)
numL =  3; % Nr. of low-spontaneous    rate (default=3)

freq2show = 1000;

%%% Stimulus generation:
fs = 100e3; % 100 kHz, sampling frequency
dt = 1/fs;
L  = [40 70 100]; % dB SPL, three stimulus levels, you can add/remove levels here
p0 = 2e-5; % Pa, reference pressure

dur = 50e-3; % 50 ms, stimulus duration
dur_click = 100e-6; % 100 us, click duration
dur_click_samples = round(dur_click*fs); % click duration in samples

t = 0:dt:dur-dt;
N_samples = length(t);
N_levels  = length(L);

insig=zeros(N_samples,N_levels); % Memory allocation. Each column is one test signal

for j=1:N_levels
    L_in_dB_peSPL = p0*10^(L(j)/20)*sqrt(2)*2; % conversion from dB to dB peSPL
    insig(10:10+dur_click_samples,j)=L_in_dB_peSPL; % one (positive) polarity click
end
%%% End Stimulus generation

%% set synaptopathy and store the results
% decide which responses you want to store: 
    % e= emission 
    % v= velocity 
    % i=ihc 
    % h=hsr 
    % m=msr 
    % l=lsr 
    % b=summed AN responses for each CF as well as CN and IC responses 
    % w=population response waves I, III, V

%% start the simulation
output = verhulst2018(insig,fs,fc_flag, ...
           'v','oae','anfH','anfM','anfL',... % model flags, specify the required output 
           'numH',numH,'numM',numM,'numL',numL); % number of AN fibres

cf = output(1).cf; % characteristic frequencies of the simulated cochlear sections
N_CFs = length(cf);

%%%% Figure 1"    CF=output.cf;
% cf is sorted in a base-to-apex (from high-to-low frequencies):
figure(1)
plot(cf), hold on
xlabel('Cochlear Channel Number [-]')
ylabel('Characteristic Frequency [Hz]')
%%%

fs_bm = output.fs_bm; %the sampling frequency of BM, OAE and IHC are higher to avoid numerical errors (see Altoe et al., 2014 JASA)
fs_an = output.fs_an; %the sampling frequency of AN, CN, IC and waves are 5 times lower
N_samples_an = size(output(1).an_summed,1);
t_an  = [0:N_samples_an-1]/fs_an; % time after the AN. By default is a decimated  version of t_bm
f = 0:fs_bm/N_samples:fs_bm-1/N_samples; % Frequency if FFT to v_BM is assessed

% Pick a channel number to plot results from. The CF corresponding to the
% channel number depends on whether you chose 'all', 'half' or 'abr':
idx_bin = find(freq2show < cf, 1,'last'); % looks for bin that is at 'freq_to_show'

% Reorganisation of the data for easier processing:
oae=zeros(N_samples,N_levels);
v=zeros(N_samples,N_levels);
ihc=zeros(N_samples,N_levels);
vrms=zeros(1,N_CFs);
ihcrms=zeros(1,N_CFs);
for j=1:N_levels
   oae(:,j)    = output(j).oae(:);
   v(:,j)      = output(j).v(:,idx_bin);
   ihc(:,j)    = output(j).ihc(:,idx_bin);
   vrms(j,:)   = rms(output(j).v);
   ihcrms(j,:) = rms(output(j).ihc);
end

%% Ear-canal pressure and Otoacoustic emissions (OAE):
% This figure is the ear-canal pressure.
% For OAE simulations:
%    - To get the reflection-source OAE, do the following simulation (see 
%      also demo_verhulst2012.m):
%          OAE_{reflections on}-OAE_{reflections off}
%    - To get the distortion-source OAE, do a normalisation using a low level
%      (linear) simulation that the reflections off

%%%% Figure 2
figure(2),
subplot(2,1,1)
plot(1000*t,oae), hold on
xlabel('Time [ms]'),ylabel('Ear Canal Pressure [Pa]'),
xlim([0 20]),
% ylim([-0.02 0.02]),
legend(num2str(L')),legend('boxoff')

subplot(2,1,2), % figure
plot(f/1000,20*log10(abs(fft(oae/p0)))); hold on
xlabel('Frequency [kHz]'),
ylabel('EC Magnitude [dB re p0]'),
xlim([0 12]),
legend(num2str(L')),
legend('boxoff')

%% v_bm and V_IHC 
figure(3),
subplot(2,2,1),
plot(1000*t,v), hold on
xlabel('Time [ms]');
ylabel('v_{bm} [m/s]');
xlim([0 30]);
title(['CF of ',num2str(round(cf(idx_bin))),' Hz']);
legend(num2str(L'));
legend('boxoff');

subplot(2,2,2),
plot(cf/1000,20*log10(vrms)), hold on
xlabel('CF [kHz]');
ylabel('rms of v_{bm} [dB re 1 m/s]');
xlim([0 14]),
% ylim([max(max(20*log10(vrms)))-100 max(max(20*log10(vrms)))+10])
title(['Excitation Pattern'])
legend(num2str(L'));
legend('boxoff');

subplot(2,2,3),
plot(1000*t,ihc), hold on
xlabel('Time [ms]');
ylabel('V_{ihc} [V]');
xlim([0 30]);
title(['CF of ',num2str(round(cf(idx_bin))),' Hz']);
legend(num2str(L'));
legend('boxoff');

subplot(2,2,4),
plot(cf/1000,ihcrms), hold on
xlabel('CF [kHz]');
ylabel('rms of V_{ihc} [V]');
xlim([0 14]);
%ylim([max(max(20*log10(ihcrms)))-100 max(max(20*log10(ihcrms)))+10])
title('Excitation Pattern');
legend(num2str(L'));
legend('boxoff');

% Reorganisation of the data for easier processing
for j = 1:N_levels
    HSR(:,j)=output(j).anfH(:,idx_bin);
    MSR(:,j)=output(j).anfM(:,idx_bin);
    LSR(:,j)=output(j).anfL(:,idx_bin);
    AN(:,j)=output(j).an_summed(:,idx_bin);
    CN(:,j)=output(j).cn(:,idx_bin);
    IC(:,j)=output(j).ic(:,idx_bin);
    W1(:,j)=output(j).w1(:);
    W3(:,j)=output(j).w3(:);
    W5(:,j)=output(j).w5(:);
    EFR(:,j)=output(j).w1+output(j).w3+output(j).w5;
end

%single unit responses
figure(4),
subplot(3,2,1),
plot(1000*t_an,HSR), hold on
title(['CF of ',num2str(round(cf(idx_bin))),' Hz'])
xlim([0 20]);
xlabel('Time [ms]');
ylabel('HSR fiber [spikes/s]')
legend(num2str(L'));
legend('boxoff')

subplot(3,2,3),
plot(1000*t_an,MSR), hold on
xlim([0 20]);
xlabel('Time [ms]');
ylabel('MSR fiber [spikes/s]')

subplot(3,2,5),
plot(1000*t_an,LSR), hold on
xlim([0 20]);
xlabel('Time [ms]');
ylabel('LSR fiber [spikes/s]');

subplot(3,2,2),
plot(1000*t_an,AN), hold on
title(['CF of ',num2str(round(cf(idx_bin))),' Hz'])
xlim([0 20]);
xlabel('Time [ms]');
ylabel('sum AN [spikes/s]');

% Spikes summed across all fibers @ 1 CF
subplot(3,2,4),
plot(1000*t_an,CN), hold on
xlim([0 20]);
xlabel('Time [ms]');
ylabel('CN [spikes/s]')

subplot(3,2,6),
plot(1000*t_an,IC), hold on
xlim([0 20]);
xlabel('Time [ms]');
ylabel('IC [spikes/s]')

% Population responses
figure(5)
subplot(4,1,1),
plot(1000*t_an,1e6*W1), hold on
title('Population Responses summed across simulated CFs');
xlim([0 20]);
xlabel('Time [ms]');
ylabel('W-1 [\muV]');
legend(num2str(L'));
legend('boxoff')

subplot(4,1,2),
plot(1000*t_an,1e6*W3), hold on
xlim([0 20]);
xlabel('Time [ms]');
ylabel('W-3 [\muV]')

subplot(4,1,3),
plot(1000*t_an,1e6*W5), hold on
xlim([0 20]);
xlabel('Time [ms]');
ylabel('W-5 [\muV]');

subplot(4,1,4),
plot(1000*t_an,1e6*EFR), hold on
xlim([0 20]);
xlabel('Time [ms]');
ylabel('EFR [\muV]');
%%%%%%%%%%%%%%%%%

amt_disp(['Showing results for frequency of ' num2str(freq2show) ' Hz']);



