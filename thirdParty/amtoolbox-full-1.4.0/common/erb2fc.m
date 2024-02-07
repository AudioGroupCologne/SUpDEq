function fc = erb2fc(erb)
%ERB2FC calculates center frequency from a given erb index
%
%   Usage: fc = erb2fc(erb)
%
%   For a given ERB index in [cam], this function calculates
%   an equivalent center frequency.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/erb2fc.php


%   #Author: The AMT Team (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

    fc = (10.^(erb/21.366)-1)./0.004368;
end


