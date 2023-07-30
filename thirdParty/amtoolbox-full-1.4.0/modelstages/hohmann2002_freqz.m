function h = hohmann2002_freqz(obj, z)
%HOHMANN2002_FREQZ  Frequency response of hohmann2002 filter or filterbank
%   Usage: h = hohmann2002_freqz(filter, z)
%          H = hohmann2002_freqz(fb, z)
%
%
%   Input parameters:
%     filter  : A filter created by HOHMANN2002_FILTER.
%     fb      : A filterbank created by HOHMANN2002.
%     z       : A vector of frequencies in z-plane for which the frequency 
%               response will be computed. z = exp(2i pi f/fs)
%
%   Output parameters:
%     h       : The complex frequency response at z. Each column represents
%               a response of a filter.
%
%     h = HOHMANN2002_FREQZ(filter, z) returns the frequency response of
%     filter created by HOHMANN2002_FILTER.
%    
%     h = HOHMANN2002_FREQZ(fb, z) returns the frequency responses of the
%     individual filters in the filterbank fb created by HOHMANN2002. The
%     responses of the individual filters are stored in columns of h. 
%
%   See also: exp_hohmann2002 demo_hohmann2002
%
%   References:
%     V. Hohmann. Frequency analysis and synthesis using a gammatone
%     filterbank. Acta Acustica united with Acoustica, 88(3):433--442, 2002.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hohmann2002_freqz.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified  
%   #Author   : Universitaet Oldenburg, tp (2002 - 2007)
%   #Author   : Piotr Majdak (2016)
%   Adapted from function gfb_*_zresponse

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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



