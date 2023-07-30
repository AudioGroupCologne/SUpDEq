function [delta,sym,asymR,asymL]=data_glasberg1990(varargin)
%DATA_GLASBERG1990 Notched-noise data for the ERB scale
%   Usage: [delta,sym,asymR,asymL]=data_glasberg1990;
%
%   DATA_GLASBERG1990 returns data from Fig.5 (left panel only) from Glasberg &
%   Moore (1990) showing notched-noise data measured both in symmetrical
%   and asymetrical conditions at a center frequency (fc) of 100 Hz.
%
%   Examples:
%   ---------
%
%   To plot Fig. 5 use :
%
%     data_glasberg1990('plot');
% 
%   References:
%     B. R. Glasberg and B. Moore. Derivation of auditory filter shapes from
%     notched-noise data. Hearing Research, 47(1-2):103--138, 1990.
%     
%     B. Moore, R. Peters, and B. Glasberg. Auditory filter shapes at low
%     center frequencies. J. Acoust. Soc. Am., 88(1):132--140, 1990.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_glasberg1990.php


%   #Author: Peter L. Soendergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% TODO: explain Data in description;

% Parse input options
definput.flags.plot = {'no_plot','plot'};
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% Plot fig.5 of Glasberg & Moore (1990)
x_delta = [0 .1 .2 .3 .4 .5];
y_sym = [87.5 85.63 84.38 81.25 78.13 73.75];
y_asymR = [83.75 80 75.63 70.63 63.75];
y_asymL = [85.63 83.75 78.13];

if flags.do_plot
    figure;
    plot(x_delta,y_sym,'*k'),hold on
    plot(x_delta(2:end),y_asymR,'>b')
    plot(x_delta(2:end-2),y_asymL,'<r')
    hold off
    axis([-.05 .55 61 94])
    xlabel('Deviation of nearer noise edge \Delta')
    ylabel('Signal level at threshold (dB SPL)')
    title('Notched-noise data from Moore et al. (1990), one subject')
    legend('Symmetric notch','Rightward asymmetric notch','Leftward asymmetric notch')
end;


