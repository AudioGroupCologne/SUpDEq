%DEMO_CULLING2004 Demo for testing culling2004.m
%
%   DEMO_CULLING2004 outputs the binaural masking level difference (BMLD) as 
%   a function of the maximum of the interaural coherence for a target phase
%   of pi.
%
%   Figure 1: Binaural masking level difference as a function of the maximum of the interaural coherence
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_culling2004.php


%   #Author : Clara Hollomey (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

fc = 1000;
coherence = 0.1:0.1:1;
phase_target = pi;
phase_int = 0;

for ii = 1: length(coherence)
  [bmld_out(ii)] = culling2004(coherence(ii),phase_target,phase_int,fc);
end

figure
plot(coherence, bmld_out)
xlabel('Maximum of IC')
ylabel('Binaural masking level difference [dB]')

