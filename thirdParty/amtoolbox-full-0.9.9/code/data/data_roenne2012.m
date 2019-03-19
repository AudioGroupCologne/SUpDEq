function [ur,fs]  = data_roenne2012(varargin)
%DATA_ROENNE2012 Unitary response
%   Usage: ur = data_roenne2012()
%
%   [ur,fs]=DATA_ROENNE2012 returns the unitary response from Roenne
%   (2012) and its sampling frequency, fs=30000.
%
%   Examples:
%   ---------
%
%   The first plot shows the unitary response in the time-domain:
%
%     [ur,fs]  = data_roenne2012;
%     plot((0:length(ur)-1)/fs,ur);
%     xlabel('Time / seconds');
%     ylabel('Amplitude');
%
%   The second plot shows the magnitude response of the unitary
%   response, normalized so the highest peak reaches 0-dB:
%
%     [ur,fs]  = data_roenne2012;
%     magresp(ur,fs,90,'fir','1');
%
%   References:
%     F. M. RÃ¸nne, T. Dau, J. Harte, and C. Elberling. Modeling auditory
%     evoked brainstem responses to transient stimuli. The Journal of the
%     Acoustical Society of America, 131(5):3903-3913, 2012. [1]http ]
%     
%     References
%     
%     1. http://scitation.aip.org/content/asa/journal/jasa/131/5/10.1121/1.3699171
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_roenne2012.php

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

% TODO: explain "unitary response"
  
s=amt_load('roenne2012','ur.mat');
ur=s.ur;

fs=30000;



