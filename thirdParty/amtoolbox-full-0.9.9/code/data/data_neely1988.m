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
%     'noplot'  Don't plot, only return data. This is the default.
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
%     Acoust. Soc. Am., 83(2):652-656, feb 1988.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_neely1988.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

L = L/100; % Stimulus level, dB SPL

for ii = 1:length(L)
  tau(ii,:) = 5 + 12.9 * 5.^(-L(ii)) * (F/1000).^(-.413);
end

