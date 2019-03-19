function data = data_boyd2012
%DATA_BOYD2012 - Data from Boyd et al. (JASA-EL, 2012)
%   Usage: data = data_boyd2012
%
%   Mean externalization scores of NH listeners extracted from top panels 
%   (1 talker condition) of Fig. 1 
%
%   Output parameters:
%     data    : structure with fields
%                 reference
%                 mix
%                 ITE with subfields BB and LP
%                 BTE with subfields BB and LP
%
%   References:
%     A. W. Boyd, W. M. Whitmer, J. J. Soraghan, and M. A. Akeroyd. Auditory
%     externalization in hearing-impaired listeners: The effect of pinna cues
%     and number of talkers. J. Acoust. Soc. Am., 131(3):EL268-EL274, 2012.
%     [1]arXiv | [2]www: ]
%     
%     References
%     
%     1. http://arxiv.org/abs/%20http://dx.doi.org/10.1121/1.3687015
%     2. http://dx.doi.org/10.1121/1.3687015
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_boyd2012.php

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

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

data.reference = 100;
data.mix = 100:-25:0;
data.ITE.BB = [100,95,58,45,33];
data.BTE.BB = [63,66,38,38,24];
data.ITE.LP = [72,75,42,26,12];
data.BTE.LP = [61,62,33,24,16];

end
