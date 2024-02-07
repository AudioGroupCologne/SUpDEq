function [bmld_out] = bmld(coherence,phase_target,phase_int,fc)
%BMLD Returns the binaural masking level difference
%
%   Usage: 
%     [bmld_out] = bmld(coherence,phase_target,phase_int,fc)
%
%
%   Input parameters:
%     coherence     : interaural coherence
%     phase_target  : phase of the target signal
%     phase_int     : phase of the interferer
%     fc            : center frequency [Hz]
%
%   BMLD returns the frequency-dependent binaural masking level
%   difference as a function of the interaural coherence and the target's
%   and interferer's phase relation
%
%   reference to be cited for this bmld formula:
%   J. F. Culling, M. L. Hawley, and R. Y. Litovsky (2005) "Erratum: The role of
%   head-induced interaural time and level differences in the speech reception
%   threshold for multiple interfering sound sources," J. Acoust. Soc. Am. 118(1), 552.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/bmld.php


%   #Author : Matthieu Lavandier

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


k = (1 + 0.25^2) * exp((2*pi*fc)^2 * 0.000105^2);
bmld_out = 10 * log10 ((k - cos(phase_target-phase_int))/(k - coherence));
if bmld_out < 0;
    bmld_out = 0;
end
%return


