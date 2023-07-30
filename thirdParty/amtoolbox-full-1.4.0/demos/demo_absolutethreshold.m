%DEMO_ABSOLUTETHRESHOLD  Standards for absolute threshold of hearing
%
%   This demos generates a simple figure that shows the behaviour of
%   the different standards for absolute thresholds of hearing.
%
%   Figure 1: Thresholds of hearing by standard
%
%      The figure shows the behaviour of the different absolute thresholds of
%      hearing.
%
%   Figure 2: High frequencies
%
%      Absolute thesholds for the ER2A and the Sennheiser HDA-200
%      earphones are provided up to 16 kHz.
%
%   See also:  absolutethreshold
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_absolutethreshold.php


%   #Author: Peter L. Soendergaard (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

amt_disp('Type "help demo_absolutethreshold" to see a description of this demo...', 'progress');

figure(1);

flow=125;
fhigh=8000;
plotpoints=30;

xrange=linspace(flow,fhigh,plotpoints);


types   = {'iso226_2003','map','er3a','er2a','hda200'};
legends = {'THR (free field)','MAP (free-field)','THR (ER-3A)','THR (ER-2A)','THR (HDA 200)'};
symbols = {'r-'         ,'k-' ,'g:'  ,'b-.','y--'  };
lines   = {2            , 1   , 2    , 2   , 2   };
ticks   = [125,250,500,1000,2000,4000,8000];

hold all;
for ii=1:numel(types)
  semiaudplot(xrange,absolutethreshold(xrange,types{ii}),...
		'tick', ticks,...
		'opts',{symbols{ii},'Linewidth',lines{ii}});
end;
legend(legends{:},'Location','North');
xlabel('Frequency (Hz)');
ylabel('Absolte threshold (dB SPL)');
xlim([2 35]);


figure(2);

flow=125;
fhigh=16000;
plotpoints=30;

xrange=linspace(flow,fhigh,plotpoints);
types   = {'er2a','hda200','iso226_2003','map','er3a','er2a','hda200'};
% legends = {'er2a','hda200'};
symbols = {'b+','y*','k' ,'r--' ,'g' ,'b:','y'};

hold all;
for ii=1:2
  semiaudplot(xrange,absolutethreshold(xrange,types{ii}),...
              'tick', ticks,...
              'opts',{symbols{ii},'Linewidth',lines{ii}});
end;
hold off;
% legend(legends{:},'Location','North');
% xlabel('Frequency (Hz)');
% ylabel('Absolte threshold (dB SPL)');
% xlim([2 41]);


% types = {'iso226_2003','map','er3a','er2a','hda200'};
% symbols = {'k' ,'r--' ,'g' ,'b:','y' };  
fc=125:125:8000; hold on; box on;  
for ii=3:numel(types),  opt={symbols{ii}, 'LineWidth', 3}; 
    semiaudplot(fc,absolutethreshold(fc,types{ii}),'opts',opt);  
end;
% legend(types); 
legend(types,'Location','North');
xlabel('Frequency (Hz)');  
ylabel('Absolute hearing threshold (dB re 20 µPa)');

