function mixer = hohmann2002_mixer(fb, delay, iterations)
%HOHMANN2002_MIXER Create new mixer object within HOHMANN2002 filterbank framework
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hohmann2002_mixer.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified  
%   #Author   : Universitaet Oldenburg, tp (2002 - 2007)
%   #Author   : Piotr Majdak (2016)
%   Adapted from function gfb_mixer_new

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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


