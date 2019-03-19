function h = hohmann2002_freqz(obj, z)
%hohmann2002_freqz  Frequency response of hohmann2002 filter or filterbank
%   Usage: h = hohmann2002_freqz(filter, z)
%          H = hohmann2002_freqz(fb, z)
%
%   h = HOHMANN2002_FREQZ(filter, z) returns the frequency response of
%   filter created by HOHMANN2002_FILTER.
%    
%   h = HOHMANN2002_FREQZ(fb, z) returns the frequency responses of the
%   individual filters in the filterbank fb created by HOHMANN2002. The
%   responses of the individual filters are stored in columns of h. 
%
%   Input parameters:
%     filter  : A filter created by HOHMANN2002_FILTER.
%     fb      : A filterbank created by HOHMANN2002.
%     z       : A vector of frequencies in z-plane for which the frequency 
%               response will be computed. z = exp(2i*pi*f/fs)
%
%   Output parameters:
%     h       : The complex frequency response at z. Each column represents
%               a response of a filter.
%
%   See also: exp_hohmann2002 demo_hohmann2002
%
%   References:
%     V. Hohmann. Frequency analysis and synthesis using a gammatone
%     filterbank. Acta Acustica united with Acoustica, 88(3):433-442, 2002.
%     
%     
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/hohmann2002_freqz.php

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

% author   : Universitaet Oldenburg, tp (Jan & Nov 2006, Jan Feb 2007)
% Adapted to AMT (PM, Jan 2016) from functions gfb_*_zresponse

if ~isfield(obj,'type'), error('Type of the object missing'); end
switch(obj.type)
  case 'gfb_Filter'
    h = (1 - obj.coefficient ./ z) .^ -obj.gamma_order * obj.normalization_factor;
  case 'gfb_analyzer'
    number_of_bands = length(obj.center_frequencies_hz);
    z = z(:);
    h = ones(length(z), number_of_bands);

    for band = 1:number_of_bands
      filter = obj.filters(band);
      h(:,band) = (1 - filter.coefficient ./ z) .^ -filter.gamma_order * filter.normalization_factor;
    end
  otherwise
    error('Unknown type of HOHMANN2002 filter object');
end


