%DEMO_LYON2011_COMPRESSIVEFUNCTIONS compressive Input-Output (I/O) function of the CARFAC model
% 
%   The I/O function of the 'cascade of asymmetric resonators with fast-acting 
%   compression' (CARFAC) modelis calculated by the RMS of the out in response 
%   tones at intensities from 10 to 100 dB SPL at the CFs of 0.5, 1, 2 and 4 kHz. 
%   The nonlinear compressive shape of the I/O function represents the cochlear 
%   mechanical compression.
%
%   Figure 1: The Input-Output (I/O) functions at a center frequency of 500 Hz
%
%   Figure 2: The Input-Output (I/O) functions at a center frequency of 1 kHz
%
%   Figure 3: The Input-Output (I/O) functions at a center frequency of 2 kHz
%
%   Figure 4: The Input-Output (I/O) functions at a center frequency of 4 kHz
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_lyon2011_compressivefunctions.php


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

%% Generate t deviation from the model assumptions 
Fs= 22050; % this is the default. You can change it by calling the CARFAC_Design function.
T= 0.1; % 100 ms length
t=0:(1/Fs):T;
fsig = [500,1000,2000,4000]; % the frequency of the tone
ch_CF=[49,39,29,16]; % the corresponding CF channels

spl=10:10:100; % intensities from 10 to 100 dB SPL
level = 20e-6.*(10.^(spl./20));
out=zeros(numel(fsig),numel(spl));

for j=1:numel(fsig)
    for i=1:length(spl)
        sig = sin(2*pi*fsig(j).*t) * level(i);
        [sig] = fade( sig,0.01,Fs ); % use 10-ms ramps to minmize the spectral splatter
        % feed the tone into the CARFC model
        [CF, decim_naps, naps, BM, ohc, agc] = lyon2011(sig', CF_struct);
        %[CF, decim_naps, naps, BM, ohc, agc] = lyon2011(CF_struct, sig');
        % calculate the RMS at the output
        out(i,j)=sqrt(2)*rms(BM(:,ch_CF(j)));
    end
end

%% Normalize the I/O function to point (10,10)
norm_out=zeros(size(out));
for i=1:length(fsig)
   out_dB=20.*log10(out./2e-5); % convert to dB SPL scale
   norm_out(:,i)=out_dB(:,i)-out_dB(1,i)+10; % Normalize
end

%% Physiological and psychoacoustic data at CF=500 Hz
% physiological chinchila data from Rhode and Cooper(1996)
data_rhode = data_lyon2011('rhode1996');

%% psychoacoustic data at C F= 500Hz, 1kHz, 2kHz, 4kHz (Lopez-Poveda et al., 2003)
data_lopez = data_lyon2011('lopezpoveda2003');

%% physiological chinchila data by Russel and Nilsen (1997)
data_russel = data_lyon2011('russel1997');

%% Illustrate the I/O functions along with the corresponding experimental data. 
figure, % CF= 500 Hz
plot(spl,norm_out(:,1));
hold on
plot(data_rhode.L_animal_500hz,data_rhode.IO_animal_norm_500hz,'xr');
hold on
plot(data_lopez.L_psych_500hz,data_lopez.IO_ex_norm_500hz,'or');
title('Input/Output function at CF = 500 Hz');
xlabel('Input [dB SPL]');
ylabel('Normalized Output [dB]');
legend('CARFAC', 'Physiological data: Rhode and Cooper, (1996)', 'Psychoacoustic data: Lopez-Poveda et al. (2003)');

figure, %CF= 1 kHz
plot(spl,norm_out(:,2)); 
hold on
plot(data_lopez.L_psych_1khz,data_lopez.IO_ex_norm_1khz,'or');
title('Input/Output function at CF = 1 kHz');
xlabel('Input [dB SPL]');
ylabel('Normalized Output [dB]');
legend('CARFAC', 'Psychoacoustic data: Lopez-Poveda et al. (2003)');

figure, % CF= 2 kHz;
plot(spl,norm_out(:,3)); 
hold on
plot(data_lopez.L_psych_2khz,data_lopez.IO_ex_norm_2khz,'or');
title('Input/Output function at CF = 2 kHz');
xlabel('Input [dB SPL]');
ylabel('Normalized Output [dB]');
legend('CARFAC', 'Psychoacoustic data: Lopez-Poveda et al. (2003)');

figure, % CF= 4 kHz;
plot(spl,norm_out(:,4));
hold on
plot(data_russel.L_animal_4khz,data_russel.IO_animal_norm_4khz,'xr');
hold on
plot(data_lopez.L_psych_4khz,data_lopez.IO_ex_norm_4khz,'or');
title('Input/Output function at CF = 4 kHz');
xlabel('Input [dB SPL]');
ylabel('Normalized Output [dB]');
legend('CARFAC', 'Physiological data: Russel and Nilsen, (1997)', 'Psychoacoustic data: Lopez-Poveda et al. (2003)');


