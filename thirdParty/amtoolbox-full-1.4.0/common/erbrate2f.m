function f = erbrate2f(erb)
%ERBRATE2F calculates the erb rate for a given frequency
%   Usage: f = erbrate2f(erb)
%
%   ERBRATE2F accepts an erb rate value as an input
%   and calculates its frequency in Hertz. It corresponds to
%   the inverse equation of Moore and Glasberg (1983) Eq. 3
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/erbrate2f.php


%   #Author: John Culling (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

f = (0.312 - (exp((erb - 43)/11.17)) * 14.675) / (exp((erb - 43)/11.17) - 1);
f = f * 1000; %convert f to Hz not kHz


