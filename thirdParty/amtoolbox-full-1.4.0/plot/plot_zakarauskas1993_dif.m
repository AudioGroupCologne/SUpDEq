function [h] = plot_zakarauskas1993_dif( d,c,pol,pos,tit )
% PLOT_ZAKARAUSKAS1993_DIF plots summed differences of estimation
%
%   Usage:
%     [h] = plot_zakarauskas1993_dif( d,c,pol,pos )
%     [h] = plot_zakarauskas1993_dif( d,c,pol,pos,tit )
%
%   Input arguments:
%     d:       summed differences of D-estimator (1st order differential)
%     c:       summed differences of C-estimator (2nd order differential)
%     pol:     source angles
%     pos:     pol-index of plotted response angle
%     tit:     string for figure title
%
%   Definition:
%     plots summed differences of estimation (Zakarauskas 1993, Fig. 9)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_zakarauskas1993_dif.php


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
if ~exist('pos','var')
    pos=23;
end
if ~exist('tit','var')
    tit='comparison of D- and C-estimator';
end

h=figure('Name','Localization model of Zakarauskas et al.(1993), cf. Fig.9','NumberTitle','off');
clf
plot(pol,d(:,pos),':');
hold on
plot(pol,c(:,pos),'r');
legend('D_n','C_n');
xlabel('Source elevation (degrees)');
ylabel('Summed difference');
title(tit);

end


