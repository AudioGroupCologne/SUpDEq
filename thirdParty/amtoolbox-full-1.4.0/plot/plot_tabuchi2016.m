function plot_tabuchi2016(gammaVal, Kave, ThreshIn)
%PLOT_TABUCHI2016 plots the lower panel from Tabuchi et al
%
%   Usage: plot_tabuchi2016(gammaVal, Kave, ThreshIn)
%
%   Input parameters:
%     gammaVal :   gamma value
%     Kave     :   K value
%     ThreshIn :   input threshold
%
%   this function plots the K value for a given gamma value and
%   input threshold
%
%   See also: tabuchi2016 exp_tabuchi2016
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_tabuchi2016.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Verified
%   #Requirements: 
%   #Author: Hisaaki Tabuchi
%   #Author: Clara Hollomey (adaptations for AMT)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

Cvec = -1.5:0.25:1; 

MaxPos = find(Cvec==-1.5);
MinPos = find(Cvec==0);
Thresh_OffFreqPre_60dB = ThreshIn(MaxPos,:,1,1) - ThreshIn(MinPos,:,1,1); % Difference of threshs of C=-1.5, and C=0
Thresh_OffFreqPre_90dB = ThreshIn(MaxPos,:,1,2) - ThreshIn(MinPos,:,1,2); 
Thresh_OnFreqPre_60dB = ThreshIn(MaxPos,:,2,1) - ThreshIn(MinPos,:,2,1); 
Thresh_OnFreqPre_90dB = ThreshIn(MaxPos,:,2,2) - ThreshIn(MinPos,:,2,2); 
xvals = [1:size(ThreshIn,2)]-1; % -1 is the adjustment because gmax starts with zero

figure
subplot(1,2,1)
%plot(xvals,Thresh_OffFreqPre_60dB,'-ok')
plot(xvals,Thresh_OffFreqPre_60dB,'m')
hold on
%plot(xvals,Thresh_OffFreqPre_90dB,'-sr')
plot(xvals,Thresh_OffFreqPre_90dB,'c','Linewidth', 2)
legend('60-dB Masker, Off-Freq', '90-dB Masker, Off-Freq')
xlabel('Gmax','FontSize',10)
ylabel('MMD; Predicted threshold difference of C=0 and C=-1.5  (dB)','FontSize',10)
%title(['GrandMean, Glasberg and Boltzmann function (beta=1, gamma=' num2str(gammaVal) ', Xshift=-20)' ', K=' num2str(Kave) ' at Gmax=34'])
xt = 0:10:70;
yt = 0:2:20;
axis([0 70 min(yt) max(yt)])
set(gca,'XTick',xt,'YTick',yt)
ylim([0 14])


subplot(1,2,2)
%plot(xvals,Thresh_OnFreqPre_60dB,'-db')
plot(xvals,Thresh_OnFreqPre_60dB,'m')
hold on
%plot(xvals,Thresh_OnFreqPre_90dB,'-^m')
plot(xvals,Thresh_OnFreqPre_90dB,'c','Linewidth', 2)
legend('60-dB Masker, On-Freq','90-dB Masker, On-Freq');
xlabel('Gmax','FontSize',10)
ylabel('MMD; Predicted threshold difference of C=0 and C=-1.5  (dB)','FontSize',10)
%title(['GrandMean, Glasberg and Boltzmann function (beta=1, gamma=' num2str(gammaVal) ', Xshift=-20)' ', K=' num2str(Kave) ' at Gmax=34'])
xt = 0:10:70;
yt = 0:2:20;
axis([0 70 min(yt) max(yt)])
set(gca,'XTick',xt,'YTick',yt)
ylim([0 14])

