%DEMO_CHEN2011 plot equal loudness contours from chen2011 to ISO standard
%
%   DEMO_CHEN2011 contrasts the equal loudness contours obtained from the 
%   model devised by Chen et al. (2011) with those from ISO 226-2003.
%
%   Figure 1: Comparison of equal loudness contours
%  
%   See also: chen2011 exp_chen2011
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_chen2011.php


%   #Author : Zhangli Chen (2011)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

OuterEarOpt = 'FreeField';

%% 1-kHz as reference
CF1k = 1000;
LdB1k = [2 10:10:90]';
for j = 1:length(LdB1k)
    Ldn1k(j) = chen2011(CF1k,LdB1k(j), OuterEarOpt);
end

%% calculating equal loudness contours
%this is the data from iso226-2003
data = amt_load('chen2011', 'data_chen2011.mat');
freq = data.freq;
spl = data.spl;

f = freq(:,1);
LdB = -10:2:140;

for i = 1:length(f)  
    for j = 1:length(LdB)
        [Ldn(j),Exitation,Cam,CF] = chen2011(f(i),LdB(j), OuterEarOpt);          
    end
    
    LdB_sameLdn(i,:) = interp1(Ldn,LdB,Ldn1k,'linear','extrap');
end

figure;
semilogx(freq(:,1),spl(:,1)-1.6,'k--',freq(:,1),LdB_sameLdn(:,1),'k',...
    freq(:,2:end),spl(:,2:end),'k--',freq(:,2:end),LdB_sameLdn(:,2:end),'k');
legend('ISO 226-2003','This paper');
xlabel('Frequency, Hz');ylabel('LdB, dB SPL');
axis([20 20000 -10 130]);
set(findobj('FontSize',10),'FontSize',12);
set(get(gca,'XLabel'),'FontSize',12);
set(get(gca,'YLabel'),'FontSize',12);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'LineWidth',1');
text(600,44,'40 phons');



