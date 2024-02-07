function delay = hohmann2002_delay(fb, delay_samples)
%HOHMANN2002_DELAY  Create new delay object within HOHMANN2002 filterbank framework
%   Usage: delay = hohmann2002_delay(fb, delay_samples)
%
%   Input parameters:
%     fb            : The filterbank structure as returned by HOHMANN2002.
%     delay_samples : The desired group delay in samples. Must be at least 1,
%                     because of the way the phase factors are computed. Larger
%                     delays lead to better signal quality.
%   Output parameters:
%     delay         : The new delay object
%
%   HOHMANN2002_DELAY(fb, delay_samples) creates a new delay object
%   that can act as the first stage of a synthesizer that
%   resynthesizes the output of the gammatone filterbank.  The
%   purpose of the delay object is to delay the output of each band by a
%   band-dependent amount of samples, so that the envelope of the impulse
%   response of the analyzer is as large as possible at the desired delay.
%   Additionally, the delay object will multiply this delayed output with a
%   band-dependent complex phase factor, so that the real part of the
%   impulse response has a local maximum at the desired delay.  Finally, the
%   delay object will output only the real part of each band.
%
%   The phase factors are approximated numerically in this constructor,
%   using a method described in Herzke & Hohmann (2007). The
%   approximation assumes parabolic behavior of the real part of the
%   impulse response in the region of the desired local maximum: The phase
%   factors are chosen so that the real parts of the impulse response in
%   the samples directly preceeding and following the desired local
%   maximum will be equal after multiplication with the pase factor.
%
%   References:
%     T. Herzke and V. Hohmann. Improved numerical methods for gammatone
%     filterbank analysis and synthesis. Acta Acustica united with Acoustica,
%     93(3):498--500, 2007.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hohmann2002_delay.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified  
%   #Author   : Universitaet Oldenburg, tp (2002 - 2007)
%   #Author   : Piotr Majdak (2016-2021): adapted to AMT; modified to time being the first dimension
%   Adapted from function gfb_*_process

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

delay.type           = 'gfb_Delay';

fb = hohmann2002_clearstate(fb);
impulse = [1; zeros(delay_samples + 2,1)];
ir = hohmann2002_process(fb, impulse).';

number_of_bands = size(ir, 1);

[~, max_indices] = max(abs(ir(:,1:(delay_samples+1))).');

delay.delays_samples = delay_samples + 1 - max_indices;

delay.memory = zeros(number_of_bands, max(delay.delays_samples));

slopes = zeros(1, number_of_bands);

for band = 1:number_of_bands
  band_max_index = max_indices(band);
  slopes(band) = (ir(band, band_max_index+1) - ir(band, band_max_index-1));
end
slopes = slopes ./ abs(slopes);
delay.phase_factors = 1i ./ slopes;


