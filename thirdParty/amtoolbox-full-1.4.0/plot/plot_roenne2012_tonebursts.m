function plot_roenne2012_tonebursts(waveVlat,click_latency)
%PLOT_ROENNE2012_TONEBURSTS plots Rønne et al. (2012) Fig. 5
%   Usage: plot_roenne2012_tonebursts(flag)
%
%   PLOT_ROENNE2012_TONEBURSTS(waveVlat,click_latency) plots the output
%   from ROENNE2012_TONEBURSTS in a similar way as Fig. 5 from the Rønne
%   et al. (2012) ABR model.
%   
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%  
%   References:
%     J. Harte, G. Pigasse, and T. Dau. Comparison of cochlear delay
%     estimates using otoacoustic emissions and auditory brainstem responses.
%     J. Acoust. Soc. Am., 126(3):1291--1301, 2009.
%     
%     S. Neely, S. Norton, M. Gorga, and J. W. Latency of auditory brain-stem
%     responses and otoacoustic emissions using tone-burst stimuli. J.
%     Acoust. Soc. Am., 83(2):652--656, feb 1988.
%     
%     F. M. Rønne, T. Dau, J. Harte, and C. Elberling. Modeling auditory
%     evoked brainstem responses to transient stimuli. The Journal of the
%     Acoustical Society of America, 131(5):3903--3913, 2012. [1]http ]
%     
%     M. S. A. Zilany and I. C. Bruce. Representation of the vowel (epsilon)
%     in normal and impaired auditory nerve fibers: Model predictions of
%     responses in cats. J. Acoust. Soc. Am., 122(1):402--417, jul 2007.
%     
%     References
%     
%     1. http://scitation.aip.org/content/asa/journal/jasa/131/5/10.1121/1.3699171
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_roenne2012_tonebursts.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Peter L. Sondergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% Center Frequencies of stimuli
CF      = [1000, 1500, 2000, 3000 ,6000, 8000];                 

%% Plot

% plot click latencies
semilogx(10e3, click_latency(1),'ko','linewidth',2);hold on;
semilogx(10e3, click_latency(2),'kv','linewidth',2);
semilogx(10e3, click_latency(3),'kd','linewidth',2);
semilogx(10e3, click_latency(4),'k*','linewidth',2);
semilogx(10e3, click_latency(5),'kx','linewidth',2);
semilogx(10e3, click_latency(6),'kp','linewidth',2);
semilogx(10e3, click_latency(7),'k^','linewidth',2);
legend('40 dB','50 dB','60 dB','70 dB','80 dB','90 dB','100 dB');
% plot toneburst latencies
set(gca,'fontsize',12);
semilogx(CF, waveVlat(:,1),'-ko','linewidth',2);
semilogx(CF, waveVlat(:,2),'-kv','linewidth',2);
semilogx(CF, waveVlat(:,3),'-kd','linewidth',2);
semilogx(CF, waveVlat(:,4),'-k*','linewidth',2);
semilogx(CF, waveVlat(:,5),'-kx','linewidth',2);
semilogx(CF, waveVlat(:,6),'-kp','linewidth',2);
semilogx(CF, waveVlat(:,7),'-k^','linewidth',2);
% Plot Neely et al (1988) data
F       = 1e3:8e3;
tau     = data_neely1988(F,40:10:100);
semilogx(F,tau','k--');
text(10e3,8, 'Click','fontsize',12);
semilogx(12e3, [5.9, 6.6, 7.6],'ok','MarkerFaceColor','k');
text(13e3,5.88, '95.2');
text(13e3,6.59, '75.2');
text(13e3,7.6, '55.2');
% Plot Harte et al (2009) data
F           = 800:10e3;
tauHarte    = 5+11.09*(F/1000).^(-0.37)/2;
semilogx(F,tauHarte,'k:','linewidth',1.5);
% Plot init
ylabel('Latency of wave V [ms]');
xlabel('Frequency of toneburst [kHz]');
set(gca,'fontsize',12);
set(gca,'XTick',CF,'XTickLabel',CF/1000);
set(gca,'YTick',1:15,'YTickLabel',1:15);
axis([800 16000 5.5 12.5]);box on;

hold off


