%DEMO_ROENNE2012 implemented for testing purposes
%
%   Simulates a click evoked ABR (c0 of the loaded file is a click). 
%   Note that the click loaded in this example starts after 15ms. 
%   The simulated wave V latency is thus also 15 ms "late"
%
%   Figure 1: ABR response
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_roenne2012.php


%   #Author: Clara Hollomey (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

stim=data_elberling2010('stim');
stim_level = 60;
flow = 100;
fhigh = 16000;

[waveVamp, waveVlat, simpot, ANout]  = roenne2012(stim.c0,30e3,stim_level);

plot_roenne2012(stim_level,waveVamp, waveVlat, simpot, ANout, 'flow',flow, 'fhigh', fhigh);



