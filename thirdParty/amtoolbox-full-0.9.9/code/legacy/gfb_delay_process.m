function [outsig, delay] = gfb_delay_process(delay, insig)
%GFB_DELAY_PROCESS  Filterbank delay processing
%   Usage: [outsig, delay] = gfb_delay_process(delay, insig)
%
%   Input parameters:
%     delay  : A gfb_delay structure created from GFB_DELAY_NEW. The delay
%              will be returned with updated delayline states as the second
%              return parameter
%     insig  : A complex matrix containing the signal to delay.  Each row
%              corresponds to a filterbank band
%
%   Output parameters:
%     outsig : A real matrix containing the delay's output
%
%   GFB_DELAY_PROCESS(delay, insig) will delay each band (row) of the
%   input data insig by a band-dependent amount of samples, then
%   multiplied with a band-dependent complex constant.  Finally, the real
%   part of this product will be returned.
%
%   See also: gfb_delay_new
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/gfb_delay_process.php

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

% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

warning('Warning: GFB_DELAY_PROCESS will be removed in a future release. Use hohmann2002_process instead. ');

[number_of_bands, number_of_samples] = size(insig);
if (number_of_bands ~= length(delay.delays_samples))
  error('input rows must match the number of bands');
end
outsig = zeros(number_of_bands, number_of_samples);
for band = [1:number_of_bands]
  if (delay.delays_samples(band) == 0)
    outsig(band,:) = ...
        real(insig(band,:) * delay.phase_factors(band));
  else
    tmp_out = [delay.memory(band,1:delay.delays_samples(band)), ...
               real(insig(band,:) * delay.phase_factors(band))];
    delay.memory(band,1:delay.delays_samples(band)) = ...
        tmp_out(number_of_samples+1:length(tmp_out));
    outsig(band,:) = tmp_out(1:number_of_samples);
  end
end


