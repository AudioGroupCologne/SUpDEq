%DEMO_EWERT2000 transfer function of the EPSM-filterbank 
%
%   DEMO_EWERT2000 outputs the transfer function of the EPSM-filterbank 
%   as presented in Ewert & Dau (2000)
%
%   Figure 1: Squared transfer functions of the filter bank
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_ewert2000.php


%   #Author : Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


fm = 50;    % Modulation frequency
l = 2;      % Length of the signal in seconds
fs = 44100; % Sampling frequency
t = 0:1/fs:l;
n = length(t);
noise = 1-2*randn(1,n);
modnoise = noise.*(1+cos(2*pi*t*fm));


[fcs, powers, TFs, freqs] = ewert2000(noise,fs);

figure
fnts = 14;
lw = 2;
% plot(freqs,10*log10(abs(TFs(1,:))),'linewidth',lw), hold on
for k = 1:7
    plot(freqs,10*log10(abs(TFs(k,:))),'linewidth',lw),hold on
end
title('Squared transfer functions of the filterbank')
xlabel('Frequency [Hz]','FontSize',fnts)
ylabel('Filter attenuation [dB]','FontSize',fnts)
set(gca,'XScale','linear','Xtick',fcs,'FontSize',fnts,'FontWeight','b');
xlim([1 79])
ylim([-20 5])
grid on


