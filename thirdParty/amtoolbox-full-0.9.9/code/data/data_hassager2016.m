function data = data_hassager2016
%DATA_HASSAGER2016 - Data from Hassager et al. (JASA, 2016)
%   Usage: data = data_hassager2016
%
%   NH data from Hassager et al. (JASA, 2016), Fig. 6,
%   representing listeners' ratings for the dir(ect-sound)
%   condition
%
%   Output parameters:
%     data    : structure with fields
%                 B ... Bandwidth Factor (ERB)
%                 angle ... source angle (deg)
%                 rating ... Externalization Rating (dim: B x angle)
%
%   References:
%     H. G. Hassager, F. Gran, and T. Dau. The role of spectral detail in the
%     binaural transfer function on perceived externalization in a
%     reverberant environment. J. Acoust. Soc. Am., 139(5):2992-3000, 2016.
%     [1]arXiv | [2]www: ]
%     
%     References
%     
%     1. http://arxiv.org/abs/http://dx.doi.org/10.1121/1.4950847
%     2. http://dx.doi.org/10.1121/1.4950847
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_hassager2016.php

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

data.B = [nan,0.316, 0.570, 1.03, 1.85, 3.33, 6.0, 10.8, 19.5 35.0, 63.1];

data.angle = [0,50];

data.rating = ...
  [ 4.79,4.73,4.70,4.67,4.43,3.65,2.81,1.97,1.60,1.49,1.30;
    4.94,4.92,4.94,4.85,4.79,4.33,3.86,3.21,2.59,2.08,1.60 ]';

end
