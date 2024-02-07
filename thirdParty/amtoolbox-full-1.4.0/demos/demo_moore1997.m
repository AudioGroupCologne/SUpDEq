%DEMO_MOORE1997
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_moore1997.php


%   #Author: Clara Hollomey (2020): for testing

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

fs = 32000; 
t = linspace(0,1,fs);
sig = sin(2*pi*1000*t).';
inSig = scaletodbspl(sig,100, 100);  

[results] = moore1997(inSig,fs);

figure
plot(results.erbN, results.eLdB)
title('Excitation pattern')
ylabel('E/E0')
xlabel('Critical bands [cams]')

