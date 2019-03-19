function [hartestim,fs]  = data_harte2009()
%DATA_HARTE2009 Tone burst stimuli from Harte et al. (2009)
%   Usage: [hartestim,fs] = data_harte2009;
%
%   [hartestim,fs]=DATA_HARTE2009 returns the tone burst stimuli from
%   Harte et al. (2009) and the sampling frequency, fs=48000.
%
%   References:
%     J. Harte, G. Pigasse, and T. Dau. Comparison of cochlear delay
%     estimates using otoacoustic emissions and auditory brainstem responses.
%     J. Acoust. Soc. Am., 126(3):1291-1301, 2009.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_harte2009.php

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

% TODO: explain stimuli in description;

hartestim = amt_load('harte2009','stim.mat');
  
fs = 48e3;
  
