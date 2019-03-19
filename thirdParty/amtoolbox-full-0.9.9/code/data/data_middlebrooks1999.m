function data = data_middlebrooks1999
%DATA_MIDDLEBROOKS1999 Statistics about non-individualized HRTFs
%   Usage: data = data_middlebrooks1999
%
%   Output parameters:
%     data.qe_own   : quadrant error rate (QE) when localizing with own
%                     HRTFs
%     data.qe_other : QE when localizing with others' HRTFs
%     data.pe_own   : local polar RMS error (PE) when localizing with own
%                     HRTFs
%     data.pe_other : PE when localizing with others' HRTFs
%     data.pb_own   : magnitude of polar bias (PB) when localizing with own
%                     HRTFs; upper-rear quadrant excluded from analysis
%     data.pb_other : PB when localizing with others' HRTFs
%
%   DATA_MIDDLEBROOKS1999 returns statistics summary from Fig. 13
%   (Middlebrooks, 1999b) showing the effect of non-individualized HRTFs.
%
%   Statistics of those parameters are stored as .mean and .quantiles*
%   representing the arithmetic mean and {0,5,25,50,75,95,100} quantiles,
%   respectively.
%
%   References:
%     J. C. Middlebrooks. Virtual localization improved by scaling
%     nonindividualized external-ear transfer functions in frequency. J.
%     Acoust. Soc. Am., 106:1493-1510, 1999.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_middlebrooks1999.php

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

% Quantiles: {0,5,25,50,75,95,100}%
  % QE
  data.qe_own.quantiles = [0,0,1,3.5,5,13,17];
  data.qe_own.mean = 4.5;
  data.qe_other.quantiles = [7.5,8,12.5,19,27.5,38,39];
  data.qe_other.mean = 21;
  % PE
  data.pe_own.quantiles = [21,23,25,27,30,34,36];
  data.pe_own.mean = 28;
  data.pe_other.quantiles = [23,33,38,42,48,54,55];
  data.pe_other.mean = 42.5;
  % EB
  data.pb_own.quantiles = [1,2.5,6,10,13,20,25.5];
  data.pb_own.mean = 10;
  data.pb_other.quantiles = [0.5,2,7,18,29,42,52.5];
  data.pb_other.mean = 19;
end
