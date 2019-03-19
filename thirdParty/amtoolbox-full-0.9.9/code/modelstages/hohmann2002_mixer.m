function mixer = hohmann2002_mixer(fb, delay, iterations)
%hohmann2002_mixer Create new mixer object within HOHMANN2002 filterbank framework
%   Usage: mixer = hohmann2002_mixer(fb, delay)
%          mixer = hohmann2002_mixer(fb, delay, iterations)
%
%   Input parameters:
%     fb         : A filterbank structure as created by HOHMANN2002
%                  The mixer created by this function can
%                  act as part of a synthesizer that resynthesizes the output
%                  from the input signal analyzed by fb. 
%     delay      : A delay structure as created by HOHMANN2002_DELAY, together
%                  with the mixer created by this function, this delay can
%                  form a synthesizer.
%     iterations : The gain factors are approximated numerically in
%                  iterations. If this parameter is omitted, then the
%                  number of iterations from fb will be used. 
%
%   HOHMANN2002_MIXER creates a mixer object with gain factors suitable
%   to calculate a weighted sum of the bands present in the output of the
%   given delay.  The gain factors are computed using a numerical optimization
%   method described in Herzke & Hohmann (2007).
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/hohmann2002_mixer.php

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

% author: Universitaet Oldenburg, tp (Jan 2002, Jan, Sep 2003, Nov 2006, Jan 2007)
% Adapted to AMT (PM, Jan 2016) from function gfb_mixer_new

if (nargin < 4)
  iterations = fb.gaincalc_iterations;
end

mixer.type = 'gfb_mixer';
center_frequencies = fb.center_frequencies_hz;
number_of_bands = length(center_frequencies);
sampling_frequency = fb.fs;

% The center frequencies in the z plain
z_c = exp(2i * pi * center_frequencies(:) / sampling_frequency);

mixer.gains = ones(number_of_bands, 1);

% compute the frequency response of each filter (col) at the center
% frequencies of all filters (row)
  pos_f = hohmann2002_freqz(fb, z_c);
  neg_f = hohmann2002_freqz(fb, conj(z_c));

% apply delay and phase correction
for band = 1:number_of_bands
  pos_f(:,band) = pos_f(:,band) * delay.phase_factors(band) .* z_c .^ -delay.delays_samples(band);
  neg_f(:,band) = neg_f(:,band) * delay.phase_factors(band) .* conj(z_c) .^ -delay.delays_samples(band);
end

% combine responses at positive and negative responses to yield
% responses for real part.
f_response = (pos_f + conj(neg_f)) / 2;

for i = 1:iterations
  % add selected spectrum of all bands with gain factors
  selected_spectrum = f_response * mixer.gains;

  % calculate better gain factors from result
  mixer.gains = mixer.gains ./ abs(selected_spectrum);
end
mixer.gains = mixer.gains.';

