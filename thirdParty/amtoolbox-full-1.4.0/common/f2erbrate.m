function rate = f2erbrate(f)
%F2ERBRATE calculates the frequency for a given erb rate
%   Usage: rate = f2erbrate(f)
%
%   erbrate2f accepts a frequency value in Hertz as an input
%   and calculates its erb rate. It corresponds to
%   the equation of Moore and Glasberg (1983) Eq. 3, with the
%   frequency being accepted in Hertz (not kHz).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/f2erbrate.php


%   #Author: John Culling (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

f = f / 1000;
rate = 11.17 * log((f+0.312)/(f+14.675)) + 43;


