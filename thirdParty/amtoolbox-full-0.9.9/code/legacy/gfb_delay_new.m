function delay = gfb_delay_new(analyzer, delay_samples)
%GFB_DELAY_NEW  Create new delay object
%   Usage: delay=gfb_delay_new(analyzer,delay_samples)
%
%   Input parameters:
%     analyzer      : The gfb_analyzer structure as returned by
%                     GFB_ANALYZER_NEW.
%     delay_samples : The desired group delay in samples. Must be at least 1,
%                     because of the way the phase factors are computed. Larger
%                     delays lead to better signal quality
%   Output parameters:
%     delay         : The new gfb_delay object
%
%   GFB_DELAY_NEW(analyzer, delay_samples) creates a new `gfb_delay
%   object` that can act as the first stage of a synthesizer that
%   resynthesizes the output of the gammatone filterbank analyzer.  The
%   purpose of the delay object is to delay the output of each band by a
%   band-dependent ammount of samples, so that the envelope of the impulse
%   response of the analyzer is as large as possible at the desired delay.
%   Additionally, the delay object will multiply this delayed output with a
%   band-dependent complex phase factor, so that the real part of the
%   impulse response has a local maximum at the desired delay.  Finally, the
%   delay object will output only the real part of each band.
%  
%   The phase factors are approximated numerically in this constructor,
%   using a method described in Herzke & Hohmann (2007).  The
%   approximation assumes parabolic behaviour of the real part of the
%   impulse response in the region of the desired local maximum: The phase
%   factors are chosen so that the real parts of the impulse response in
%   the samples directly preceeding and following the desired local
%   maximum will be equal after multiplication with the pase factor.
%
%   References:
%     T. Herzke and V. Hohmann. Improved numerical methods for gammatone
%     filterbank analysis and synthesis. Acta Acustica united with Acoustica,
%     93(3):498-500, 2007.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/gfb_delay_new.php

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
% date     : Jan 2002; Nov 2003; Mar Jun Nov 2006; Jan Feb 2007

warning('Warning: GFB_DELAY_NEW will be removed in a future release. Use hohmann2002_delay instead. ');

delay.type           = 'gfb_Delay';

  analyzer             = hohmann2002_clearstate(analyzer);
  impulse              = zeros(1, delay_samples + 2);
  impulse(1)           = 1;

    impulse_response = ...
      gfb_analyzer_process(analyzer, impulse);

number_of_bands      = size(impulse_response, 1);

[dummy, max_indices] = max(abs(impulse_response(:,1:(delay_samples+1))).');

delay.delays_samples = delay_samples + 1 - max_indices;

delay.memory         = zeros(number_of_bands, max(delay.delays_samples));

slopes = zeros(1, number_of_bands);
for band = [1:number_of_bands]
  band_max_index = max_indices(band);
  slopes(band) = (impulse_response(band, band_max_index+1) - ...
                  impulse_response(band, band_max_index-1));
end
slopes = slopes ./ abs(slopes);
delay.phase_factors = 1i ./ slopes;


