function [weightings] = siiweightings(fc)
%SIIWEIGHTINGS  Compute the SII weightings
%   Usage: [weightings] = siiweightings(fc)
%
%    SIIWEIGHTINGS(fc) computes the SII-weighting for the centre
%   frequencies given in fc.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/siiweightings.php

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
  
weightings = zeros(length(fc),1);
bands = [0, 100, 200, 300, 400, 4400, 5300, 6400, 7700, 9500].';
weights = [0, 0.0103, 0.0261, 0.0419, 0.0577, 0.0460, 0.0343, 0.0226, 0.0110, 0].';
for n = 1:length(fc)
  if fc(n) >= 9500
    weightings(n) = 0;
  else
    ii = find(bands > fc(n));
    weightings(n) = weights((ii(1) - 1));
  end
end
weightings = weightings ./ sum(weightings);

