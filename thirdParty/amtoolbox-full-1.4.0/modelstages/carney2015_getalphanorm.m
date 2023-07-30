function [B,A] = carney2015_getalphanorm(tau, fs, t)
%CARNEY2015_GETALPHANORM  Returns filter coefficients for a normalized alpha function
%
%   Returns z-transform coefficients for the function:
%
%       y(t) = t*e^(-t/tau);
%
%   The resulting coefficents can then be used in filter().  This version
%   normalizes the alpha function so that the area from 0 to t is equal
%   to 1.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/carney2015_getalphanorm.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Authors: University of Rochester (UR EAR) team
%   #Author: Mike Anzalone (2004)
%   #Authors: Clara Hollomey (2020): integration in the AMT
%   #Authors: Piotr Majdak (2021): integration for the AMT 1.0
%   #Authors: Alejandro Osses (2021): extensions for the AMT 1.1

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


a = exp(-1/(fs*tau));
% norm = zeros(length(t));
norm = 1 ./(tau^2 .* (exp(-t/tau) .* (-t/tau-1) + 1));

B = [0 a];
A = [1 -2*a a^2] * fs * 1 ./norm;



