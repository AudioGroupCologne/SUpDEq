%DEMO_VIEMEISTER1979 test demo
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_viemeister1979.php


%   #Author: Clara Hollomey (2020): for testing

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

fs = 44100;
duration = 1;
t = 1/fs:1/fs:duration;

carrier = sind(6000 * t);
carrier = randn(1,44100);

faxis = [0 2 4 8 16 32 64 125 250 500 1000 2000 4000];

for ii = 1:length (faxis)
  insig = sind(faxis(ii) * t).* carrier;
  modelOut(ii) = viemeister1979(insig,fs);
end  

plot(faxis(2:end), 20*log10(modelOut(2:end))/2)
%m...ratio of peak value to dc
xlim([1 1000])
xticks(faxis)

set (gca (), "ydir", "reverse")
xlabel('Modulation frequency')
ylabel('Modulation index')
grid on

