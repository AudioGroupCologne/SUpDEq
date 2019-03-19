function data = data_sabin2005
%DATA_SABIN2005 - Data retrieved from Figures 8-10 of Sabin et al. (2005)
%   Usage: data = data_sabin2005
%
%   Output parameters:
%     data.gain  : gain of regression line. Field f*
%                  for frontal targets, field r for rear targets. Field
%                  m for across-listener means, field sd for standard
%                  deviation across listeners
%     data.pqv   : percent quasi-veridical responses (within +-45 deg 
%                  relative to regression line. Fields as for gain.
%     data.var   : standard deviation of quasi-veridical responses from 
%                  regression line . Fields as for gain.
%     data.SL    : sensation level in dB (re audibility threshold)
%
%   DATA_SABIN2005 returns percentage of quasi-veridical responses
%   (audible and within 45 deg of the regression line) as a function of
%   senstation level averaged across all listeners.
%
%   References sabin2005nearthreshold
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_sabin2005.php

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

% AUTHOR: Robert Baumgartner

data.SL = [0:5:20,20:10:60];

% Figure 8
data.gain.f.m =    [-0.08,0.44,0.41,0.61,0.67,0.73,0.81,0.80,0.79,0.86];
data.gain.f.sd =   [0.12,0.13,0.12,0.10,0.08,0.10,0.06,0.04,0.05,0.03];
data.gain.r.m =    [nan,0.59,0.85,0.97,0.74,0.63,0.88,0.76,0.84,0.73];
data.gain.r.sd =   [nan,nan,0.17,0.13,0.14,0.10,0.13,0.05,0.10,0.09];

% Figure 9
data.pqv.f.m =     [21.4,35.0,71.7,81.2,84.5,81.4,89.0,88.3,90.4,87.4];
data.pqv.f.sd =    [11.9,23.1,13.5,07.5,06.0,10.1,09.6,09.2,07.2,06.0];
data.pqv.r.m =     [22.9,25.8,26.2,33.4,57.7,60.4,75.0,74.8,77.1,78.8];
data.pqv.r.sd =    [04.9,14.4,20.9,19.4,14.6,15.2,16.3,19.0,16.3,20.1];

% Figure 10
data.var.f.m =    [13.6,14.0,17.0,16.3,14.4,14.4,12.5,11.0,12.3,12.9];
data.var.f.sd =   [2.85,2.21,1.00,0.43,1.21,1.21,0.57,0.85,1.07,1.14];
data.var.r.m =    [nan,nan,32.6,23.9,20.3,17.3,17.3,16.0,15.9,16.3];
data.var.r.sd =   [nan,nan,7.41,3.99,2.49,1.57,1.50,0.85,1.35,1.28];

end
