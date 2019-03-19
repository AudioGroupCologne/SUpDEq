function synth = hohmann2002_synth(fb, desired_delay)
%hohmann2002_synth  Create new synthesizer object within HOHMANN2002 filterbank framework
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
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/hohmann2002_synth.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% author: Universitaet Oldenburg, tp (Jan 2002, Jan, Sep 2003, Nov 2006, Jan 2007)
% Adapted to AMT (PM, Jan 2016) from function gfb_synthesizer_new

synth.type = 'gfb_Synthesizer';
desired_delay_in_samples = round(desired_delay * fb.fs);
if (desired_delay_in_samples < 1)
    error('delay must be at least 1 sample');
end

synth.delay = hohmann2002_delay(fb, desired_delay_in_samples);
synth.mixer = hohmann2002_mixer(fb, synth.delay);

