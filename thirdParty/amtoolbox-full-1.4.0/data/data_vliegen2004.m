function data = data_vliegen2004
%DATA_VLIEGEN2004 experimental results from vliegen2004level for stimulus duration of 3 msec
%   Usage: data = data_vliegen2004
%
%   Output parameters:
%     data          : struct containing the data
%
%   DATA_VLIEGEN2004 returns listeners' polar angle gains obtained in 
%   free-field localization experiments with 3-ms long stimuli Data was retrieved from  
%   TAB.I, FIG.6, and FIG.7 in Vliegen & van Opstal (2004) Six listeners were tested
%
%   The data struct comprises the following fields:
%
%     'id'        listener identification
%     'SL'        probed sensation level
%     'pgf'       polar angle gain for the front
%     'var'       polar angle gain for the front
%
%   References:
%     J. Vliegen and A. J. Van Opstal. The influence of duration and level on
%     human sound localization. J. Acoust. Soc. Am., 115(1705), 2004.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_vliegen2004.php


%   #Author: Robert Baumgartner

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%%
SPL = [33    43    53    58    63    68    73]; % data from session 2 and 3

data(1).id = 'JV';
data(1).SPL = SPL;
data(1).SL2SPL = 23; % data from Tab. I
data(1).SL =  SPL - data(1).SL2SPL;
data(1).var = [20.4545   15.9091    7.7273    7.2727   11.8182    7.2727   10.0000];
data(1).pgf = [0.3200    0.6800    0.7200    0.7400    0.6400    0.5600    0.5200];

data(2).id = 'MW';
data(2).SPL = SPL;
data(2).SL2SPL = 30;
data(2).SL =  SPL - data(2).SL2SPL;
data(2).var = [13.6364    9.5455    8.6364   15.4545   10.9091   12.7273   16.3636];
data(2).pgf = [0.3000    0.4800    0.4400    0.7200    0.4400    0.5600    0.4200];

data(3).id = 'HV';
data(3).SPL = SPL;
data(3).SL2SPL = 31;
data(3).SL =  SPL - data(3).SL2SPL;
data(3).var = [14.5455   12.7273   13.6364   10.9091   11.8182   12.7273   15.4545];
data(3).pgf = [0.5600    0.4800    0.5400    0.7200    0.7200    0.6400    0.6000];

data(4).id = 'MZ';
data(4).SPL = SPL;
data(4).SL2SPL = 27;
data(4).SL =  SPL - data(4).SL2SPL;
data(4).var = [15.9091   13.1818   13.6364   14.0909    9.5455    9.0909    9.0909];
data(4).pgf = [0.1400    0.3600    0.5600    0.6200    0.5200    0.5800    0.5200];

data(5).id = 'WV';
data(5).SPL = SPL;
data(5).SL2SPL = 25;
data(5).SL =  SPL - data(5).SL2SPL;
data(5).var = [14.5455    8.1818    5.4545    7.2727    7.7273   11.3636    8.1818];
data(5).pgf = [0.3600    0.5200    0.5200    0.6000    0.4200    0.4400    0.4000];

data(6).id = 'FW';
data(6).SPL = SPL;
data(6).SL2SPL = 31;
data(6).SL =  SPL - data(6).SL2SPL;
data(6).var = [14.5455   12.7273   10.9091   13.6364   10.4545   12.2727   13.1818];
data(6).pgf = [0.2000    0.3400    0.3200    0.5400    0.2600    0.3600    0.1800];

%% Plot
% 
% figure
% for ii = 1:length(data)
%   subplot(3,2,ii)
%   plot(data(ii).SL,data(ii).pgf,'ko-')
%   hold on
%   axis([1 69 0 1.2])
%   axis square
%   text(55,1,data(ii).id)
%   xlabel('SL (dB)')
%   ylabel('Elevation gain')
% end
% 
% figure
% for ii = 1:length(data)
%   subplot(3,2,ii)
%   plot(data(ii).SL,data(ii).var,'ko-')
%   hold on
%   axis([1 52 0 25])
%   axis square
%   text(40,20,data(ii).id)
%   xlabel('SL (dB)')
%   ylabel('Response variability (deg)')
% end

end


