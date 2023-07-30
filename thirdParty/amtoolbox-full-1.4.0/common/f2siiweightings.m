function [weightings] = f2siiweightings(fc)
%F2SIIWEIGHTINGS  Compute the SII weightings
%   Usage: [weightings] = siiweightings(fc)
%
%   siiweightings(fc) derives the SII-weighting for the centre
%   frequencies given in fc.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/f2siiweightings.php


%   #Author: The AMT Team (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
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


