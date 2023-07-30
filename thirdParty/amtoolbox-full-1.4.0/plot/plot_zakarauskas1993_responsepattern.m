function [h] = plot_zakarauskas1993_responsepattern( errd,errc,pol,tit )
% PLOT_ZAKARAUSKAS1993_RESPONSEPATTERN plots response pattern estimation
%
%   Usage:       
%     plot_zakarauskas1993_responsepattern( errd,errc,pol )
%     plot_zakarauskas1993_responsepattern( errd,errc,pol,tit )
%
%   Input arguments:
%     errd:    error of D-estimator (1st order differential)
%     errc:    error of C-estimator (2nd order differential)
%     pol:     source angles
%     tit:     string for figure title
%
%   Description:
%     plots response pattern of estimation according to
%     Zakarauskas et al. (1993), Fig.8
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_zakarauskas1993_responsepattern.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Robert Baumgartner (2010), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%default settings
if ~exist('tit','var')
    tit='response patterns of D- and C-estimator';
end

h=figure('Name','Localization model of Zakarauskas et al.(1993), cf. Fig.8','NumberTitle','off');
clf
plot(pol,pol+errd,'b+:');
hold on
plot(pol,pol+errc,'rs-');
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
axis square
legend('D_n','C_n','Location','SouthEast');
xlabel('Source elevation (degrees)');
ylabel('Located elevation (degrees)');
title(tit);

end


