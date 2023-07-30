%DEMO_LYON2011 excitation patterns of the CARFAC model
%
%   This script generates a tone at 30 dB SPL and feeds it to the cascade of 
%   asymmetric resonators with fast-acting compression (CARFAC) model. The 
%   excitation pattern is formed by taking the RMS of the output channels.
%
%   Figure 1: The response at channel 16 to a 4-kHz tone at 30 dB SPL
%
%   Figure 2: The excitation pattern in response to a 4-kHz tone at 30 dB SPL
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_lyon2011.php


%   #Author: Amin Saremi, PhD. (amin.saremi@uni-oldenburg.de)
%   #Author: Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

CF = lyon2011_design; % The defauls parameter designs are here
CF_struct = lyon2011_init(CF); % This initializes all the states for the CARFACT

% The following determines the lowest and highest frequencies and the
% position-frequency relation of the channels
Flow=CF_struct.pole_freqs(end);
Fhigh= CF_struct.pole_freqs(1);
N=CF_struct.n_ch;
[ Pos,CF_CARFAC ] = f2bmdistance( Flow,Fhigh,N );


%% Create the stimuli: a 100ms-long tone at 4 kHz and 30 dB SPL
Fs= 22050; % this is the default. You can change it by calling the CarfacDesign function.
T= 0.1; % 100 ms length
t=0:(1/Fs):T;
fsig = 4000; 
SPL=30; 
level = 20e-6.* SPL.^(10/20); % 30 dB SPL

sig = sin(2*pi*fsig.*t) * level;
t_fade=0.01;
[ sig ] = fade( sig,t_fade,Fs ); % ramp the signal for 10 ms
% to minimize the spectral splatter effect.

% Now feed the signal into the model.
%[CF, decim_naps, naps, BM, ohc, agc] = lyonC2011(CF_struct,sig',0,1);
[CF, decim_naps, naps, BM, ohc, agc] = lyon2011(sig', CF_struct);

% Now plot the CARFAC output at the channel corresponding to CF of 4 kHz
if fsig==4000
    ch_CF=16;
end
figure
plot(t,BM(:,ch_CF));
title('The CARFAC response at channel 16 (CF=4kHz) to a 4-kHz tone');
xlabel('Time [s]');
ylabel('Amplitude');
%% Illustare the Excitation Pattern (pre-assumption: The channels are not distorted).
Fr=zeros(1,CF_struct.n_ch);
% calculate the RMS of the output (steady state) for all channels
for i=1: CF_struct.n_ch
    Fr(i)=sqrt(2).*rms(BM(floor(length(BM(:,i))/2):end,i));
end
norm_FR=zeros(size(Fr));
norm_Fr=Fr./max(Fr);
figure
plot(Pos, db(norm_Fr));
title('The Excitation Pattern in response to 4-kHz tone at 30 dB SPL');
xlabel('Cochlear location [0=Base; 1=Apex]')
ylabel('Normalized RMS energy');

%% physiological data (Ren,2002)%C.H. I don't think this works as intended
%data_ren = data_lyon2011('ren2002');

%hold on,
%plot(data_ren.P_ex,data_ren.Ex_30dB_norm,'xr');
%hold on
%plot(data_ren.P_ex,data_ren.Ex_30dB_norm,'r');
%legend('CARFAC','Physiological Data: Ren(2002)');




