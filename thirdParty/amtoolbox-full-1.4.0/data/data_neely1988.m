function tau  = data_neely1988(F,L,varargin)
%DATA_NEELY1988 ABR wave V data as functon of level and sweeping rate
%   Usage: tau = data_neely1988(flag)
%
%   Input parameters:
%     F    : Centre frequencies of stimulus
%     L    : Levels at centre frequencies
%
%   Output parameters:
%     tau  :  Wave V latency
%
%   DATA_NEELY1988(F,L) returns data points based on equation 1 from
%   Neely et al. (1988) where F is the centre frequencies and L are
%   the associated levels. 
%
%   The flag may be one of:
%
%     'no_plot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%
%   Examples:
%   ---------
%
%   Figure 2 in Neely et al. (1988) can be reproduced using:
%
%     F=[250 500 1000 2000 5000 8000];
%     tau=data_neely1988(F,[40 60 80 100]);
%     loglog(F,tau','k-');
%     xlabel('CF');
%     ylabel('Latency [ms]')
%
%   References:
%     S. Neely, S. Norton, M. Gorga, and J. W. Latency of auditory brain-stem
%     responses and otoacoustic emissions using tone-burst stimuli. J.
%     Acoust. Soc. Am., 83(2):652--656, feb 1988.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_neely1988.php


%   #Author: Peter L. Soendergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

L = L/100; % Stimulus level, dB SPL

for ii = 1:length(L)
  tau(ii,:) = 5 + 12.9 * 5.^(-L(ii)) * (F/1000).^(-.413);
end


