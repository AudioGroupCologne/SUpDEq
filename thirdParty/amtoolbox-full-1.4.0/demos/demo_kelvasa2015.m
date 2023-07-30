DEMO_KELVASA2015

%   #Author: Clara Hollomey (2020): for testing purposes 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_kelvasa2015.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


[results] = kelvasa2015(insig,fs,display_level);

display_level = 'debug'; % set to 'debug' to see more information, set to 'no_debug' to have less mess on your display
fc = 1000;
fs = 44100;
itd = 0.001;
ild = 1;

insig = sig_itdildsin(fc,itd,ild,fs);

[results] = kelvasa2015(insig,fs,display_level);

