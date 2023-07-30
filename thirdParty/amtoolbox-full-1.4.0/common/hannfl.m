function h = hannfl(len,h1len,h2len);
%HANNFL plots a hanning window
%
%   Input parameters:
%     len   : window length [samples]
%     h1len : lenght [samples] from left minimum to window maximum
%     h2len : lenght [samples] from window maximum to right minimum
%
%   HANNFL yields an asymmetric hanning window. Its symmetry
%   can be controlled by the ration between h1len and h2len.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/hannfl.php


%   #Author: The AMT Team (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if len > 0
h = ones(len,1);

switch h1len
case 0
otherwise
   h(1:h1len)=(1-cos(pi/h1len*[0:h1len-1]))/2;
end

switch h2len
case 0
otherwise
   h(end-h2len+1:end)=(1+cos(pi/h2len*[1:h2len]))/2;
end

else
end


