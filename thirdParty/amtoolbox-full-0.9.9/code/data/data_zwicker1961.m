function data = data_zwicker1961(varargin)
%DATA_ZWICKER1961  Data for the Bark scale
%   Usage: data = data_zwicker1961;
%
%   DATA_ZWICKER1961 returns the data points defining the notion of
%   critical bands. The output data consists of the frequency limit of
%   each band in Hz.
%
%   To get the bandwidth of each channel, simply use:
%
%     bw = diff(data_zwicker1961);
%
%   The first entry has been modified from the original paper: It was
%   originally 20 Hz, but in Zwicker and Fastl 1999 this was changed to 0
%   Hz.
%
%   The data can be plotted using :
%
%     data_zwicker1961('plot');
%  
%   References:
%     E. Zwicker. Subdivision of the audible frequency range into critical
%     bands (frequenzgruppen). J. Acoust. Soc. Am., 33(2):248-248, 1961.
%     [1]http ]
%     
%     E. Zwicker and H. Fastl. Psychoacoustics: Facts and models, volume 254.
%     Springer Berlin, 1999.
%     
%     References
%     
%     1. http://link.aip.org/link/?JAS/33/248/1
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_zwicker1961.php

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

% TODO: describe Data in description
  
definput.flags.plot = {'noplot','plot'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

data = [0,100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,...
        2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500].';

if flags.do_plot
  dl=length(data);
  semilogx(data(2:end),1:dl-1);
  xlim([100, 16000]);
  xlabel('Frequency in Hz.');
  ylabel('Critical-band function (Tonheit) in Bark.');
end;
