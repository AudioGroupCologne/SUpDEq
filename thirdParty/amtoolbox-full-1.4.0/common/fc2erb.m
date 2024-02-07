function erb = fc2erb(fc)
%FC2ERB calculates the ERB index from a given frequency
%
%   Usage: erb = fc2erb(fc)
%
%   FC2ERB accepts a vector of center frequencies
%   fc as an input and converts them to their index
%   on the ERB scale in [cam].
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/fc2erb.php


%   #Author: The AMT Team (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

    erb = 21.366*log10(0.004368.*fc + 1);
end


