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
%     notched-noise data. Hearing Research, 47(1-2):103-138, 1990.
%     
%     B. Moore, R. Peters, and B. Glasberg. Auditory filter shapes at low
%     center frequencies. J. Acoust. Soc. Am., 88(1):132-140, 1990.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_glasberg1990.php

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


% TODO: explain Data in description;

% Parse input options
definput.flags.plot = {'noplot','plot'};
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

