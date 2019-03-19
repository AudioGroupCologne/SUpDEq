function [acrossSubdSRTs,std_acrossSubSRTs]  = data_joergensen2011(varargin)
%DATA_JOERGENSEN2011 Experimental data from Joergensen and Dau (2011)
%   Usage: [acrossSubdSRTs,std_acrossSubSRTs]  = data_joergensen2011(flag)
%
%   DATA_JOERGENSEN2011 returns measured change of SRT tested by Joergensen 
%   and Dau (2011). Average and standard deviation across listeners.
%
%   The flag may be one of:
%
%     'fig5' : The measured change in SRT (open squares), averaged across 6 
%              normal-hearing listeners, as a function of the reverberation 
%              time, T30 (0s, 0.4s, 0.7s, 1.3s, or 2.3s). 
%              The mean SRT in the reference condition was -3 dB. 
%              Model predictions are indicated by the filled squares. The 
%              linear correlation coefficient (q) and RMSE is indicated in 
%              the upper left corner.
%     'fig6' : DSRT (left ordinate) as a function of the over-subtraction 
%              factor a (0, 0.5, 1, 2, 4, or 8)
%              for 4 normal-hearing listeners (open squares) and 
%              sEPSM predictions (filled squares). The right ordinate (with 
%              a reversed scale) shows the corresponding sSTI values as filled 
%              gray circles. These values are, however, not converted to DSRT 
%              values since these would be outside the left ordinate scale.
%
%   References:
%     S. Joergensen and T. Dau. Predicting speech intelligibility based on
%     the signal-to-noise envelope power ratio after modulation-frequency
%     selective processing. J. Acoust. Soc. Am., 130(3):1475-1487, 2011.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_joergensen2011.php

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

% Last modified: 2016/01/15 by Robert Baumgartner

% Define input flags
definput.flags.type = {'fig5','fig6'};

[flags,kv]  = ltfatarghelper({}, definput,varargin);

if flags.do_fig5
  % The following data was retrieved from FIG. 5
  acrossSubdSRTs =    [0.0100    3.6300    4.6700    6.8700    8.2900];
  std_acrossSubSRTs = [0.8500    0.4800    0.4800    1.1600    1.9700];
end

if flags.do_fig6
  acrossSubdSRTs =    [0, 1.3750, 1.3750, 1.6250, 1.8750, 2.7083];    
  std_acrossSubSRTs = [0.4513, 0.7500, 0.2500, 0.4590, 0.3436, 0.5159];
end
