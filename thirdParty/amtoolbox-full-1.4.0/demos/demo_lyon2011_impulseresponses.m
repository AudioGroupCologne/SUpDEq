%DEMO_LYON2011_IMPULSERESPONSES impulse responses of the CARFAC model
%
%   This script produces an 80-us impulse (click) feeds it into the the 'cascade 
%   of asymmetric resonators with fast-acting compression' (CARFAC) model.
%   The model's output ('impulse response') is converted to the frequency domain
%   by taking the FFT (yielding the 'frequency response'). The ERB
%   bandwidth and the quality factor (QERB) is calcualted and compared to the
%   experimental data.
%   The impulses responses to clicks at various intensities are also assessed
%   and the QERBs are compared with the respective experimental data.
%
%   Figure 1: Impulse responses at a center frequency of 4 kHz
%
%   Figure 2: Equivalent rectangular bandwidth at low intensities
%
%   Figure 3: The filter bandwidth as a function of intensity
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_lyon2011_impulseresponses.php


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

%% impulse responses%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=[500,1000,2000,4000]; % the CFs at which the impulse responses are studied
ch_CF=[49,39,29,16]; % the channels corresponding to above CFs.

%% click response
Fs= 22050; % this is the default. You can change it by calling the CARFAC_Design function.
T_clk=50e-3; % the length of the entire clock = 50 ms
t_clk=0:1/Fs:T_clk;
% create the click: 80 us of one and the rest is zero.
clk=zeros(1,length(t_clk));
clk(round(Fs/1000)+1:1+round(Fs/1000)+80*round(Fs/1e6))=1;      % Warning: integer operands required (colon) -> schneidet nach dem Komma ab
% the intensity from 10 to 100 dB SPL
spl=10:10:100;level = 20e-6.*(10.^(spl/20));

clk_sig=zeros(length(clk)); %zero pad 
% create a 3 dimensional matrix containing the responses at all channels at all intensities.
BM_clk=zeros(length(t_clk),71,length(spl)); 
en=zeros(71,length(spl));
phase2_ch_CF=zeros(length(t_clk), length(spl)); % The phase of the click response
BW=zeros(numel(spl),numel(f)); % the bandwidth of the response at each intensity
QREB=zeros(numel(spl),numel(f));
group_delay=zeros(size(phase2_ch_CF));

for j=1:numel(f)
    for i=1:length(spl)
        clk_sig=(level(i)*2).*clk; % dB SPL to equivalent peak-to-peak SPL for click (dB pe SPL)
        % run the model and save the results in the 'BM_clk' matrix
        [CF, decim_naps, naps, BM_clk(:,:,i), ohc, agc] = lyon2011(clk_sig', CF_struct);
        % send the BM_clk at the corresponding CF to the function 'erbest' to
        % estimate the ERB bandwidth and the quality factor (QERB).
        [BW(i,j), QERB(i,j)]= erbest(BM_clk(:,ch_CF(j),i), f(j), Fs);
        % find the channel that has the maximum energy
        for k=1:71
            en(k,i)=rms(BM_clk(:,k,i));
        end
        spectrum_clk_ch_CF=fft(BM_clk(:,ch_CF(j),i)); % take the FFT at the CF channel
        phase2_ch_CF(:,i)=unwrap(angle(spectrum_clk_ch_CF)); % unwarp the phase
        % calculate the group delay
        for k=1:length(phase2_ch_CF(:,i))-1
            group_delay(k,i)=-1*(phase2_ch_CF(k+1,i)-phase2_ch_CF(k,i))/((2*pi)/T_clk);
        end   
    end
end

[max_en,ch_max_clk]=max(en(:,3));%find the best channel for intensity at spl(3).

%% Impulse responses at CF= 4 kHz
figure
subplot(3,1,1);plot(t_clk,BM_clk(:,ch_max_clk,6))
axis([0 0.01 -2 2])
legend('response to 60-dB pe SPL click')
title('Impulse Responses at CF=4 kHz to Clicks at 20, 40 and 60 dB pe SPL')
subplot(3,1,2)
plot(t_clk,BM_clk(:,ch_max_clk,4))
axis([0 0.01 -0.5,0.5])
legend('response to 40-dB pe SPL click')
subplot(3,1,3)
plot(t_clk,BM_clk(:,ch_max_clk,2))
axis([0 0.01 -0.1, 0.1])
legend('response to 20-dB pe SPL click')
xlabel(' time [s]')
%% QERB at low intensities (20 dB pe SPL)
figure
semilogy(f,QERB(2,:),'sb')
% Experimentally- derived equation of Glasberg and Moore (1990)
% the QERBs at 0.5, 1, 2 and 4 kHz at low intensities according to Glasberg and Moore (1990).
data_glasberg = data_lyon2011('glasberg1990');

hold on
semilogy(f,data_glasberg.QERB_exp,'xr')
hold on
semilogy(f,QERB(2,:),'b')
hold on
semilogy(f,data_glasberg.QERB_exp,'r')
axis([500,4000,1,10])
title('Q_{ERB} at low intensities (20 dB pe SPL)')
xlabel('CF [Hz]')
legend('CARFAC model','Experimental data: Glasberg and Moore (1990)')
%%
figure
semilogy(spl,QERB(:,4))
% Physiological data: deBoer and Nuttal (2000)
data_deboer = data_lyon2011('deboer2000');

hold on
semilogy(data_deboer.intensity_deBoer,data_deboer.QERB_deBoer, 'xr')
axis([10 100 1 100])
legend(' CARFAC model','Physiological Data: deBoer and Nuttal (2000)')
title('Q_{ERB} as a function of intensity at CF= 4 kHz')
xlabel('Click Intensity [dB pe SPL]')
ylabel('Q_{ERB}')


