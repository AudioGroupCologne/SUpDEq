function synth = hohmann2002_synth(fb, desired_delay)
%HOHMANN2002_SYNTH  Create new synthesizer object within HOHMANN2002 filterbank framework
%   Usage: synth = hohmann2002_synth(fb, desired_delay)
%
%   Input parameters:
%     fb            : The filterbank structure as returned by HOHMANN2002.
%     desired_delay : the desired group delay of the total analysis-synthesis
%                     system (in seconds). Minimum delay is 1 sample. 
%                     Greater delays result in better output signal quality. 
%
%   Output parameters:
%     synth : the constructed synthesizer structure
%
%   HOHMANN2002_SYNTH creates a new synthesizer object for the
%   reconstruction of signals analyzed by fb within the HOHMANN2002
%   filterbank framework.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hohmann2002_synth.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified  
%   #Author   : Universitaet Oldenburg, tp (2002 - 2007)
%   #Author   : Piotr Majdak (2016)
%   Adapted from function gfb_synthesizer_new

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

synth.type = 'gfb_Synthesizer';
desired_delay_in_samples = round(desired_delay * fb.fs);
if (desired_delay_in_samples < 1)
    error('delay must be at least 1 sample');
end

synth.delay = hohmann2002_delay(fb, desired_delay_in_samples);
synth.mixer = hohmann2002_mixer(fb, synth.delay);


