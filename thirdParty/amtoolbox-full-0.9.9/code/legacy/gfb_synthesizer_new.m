function synthesizer = gfb_synthesizer_new(analyzer, desired_delay_in_seconds)
%GFB_SYNTHESIZER_NEW  Create new synthesizer
%   Usage: synthesizer = gfb_synthesizer_new(analyzer, desired_delay_in_seconds)
%
%   Input parameters:
%     analyzer  :    an analyzer struct as returned by gfb_analyzer_new
%     desired_delay_in_seconds :
%                    the desired group delay of the analysis-synthesis
%                    system in seconds.  Greater delays result in better
%                    output signal quality.  Minimum delay is (1 /
%                    analyzer.fs)
%   Output parameters:
%     synthesizer : the constructed gfb_Synthesizer structure
%
%   GFB_SYNTHESIZER_NEW creates a new synthesizer object that fits to the
%   given analyzer.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/gfb_synthesizer_new.php

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

% author   : tp
% date     : Jan 2002, Nov 2003, Nov 2006, Jan 2007

warning('Warning: GFB_SYNHTESIZER_NEW will be removed in a future release. Use hohmann2002_synth instead. ');

synthesizer.type         = 'gfb_Synthesizer';
desired_delay_in_samples = round(desired_delay_in_seconds * ...
			         analyzer.fs);
if (desired_delay_in_samples < 1)
    error('delay must be at least 1/analyzer.fs');
end

synthesizer.delay = gfb_delay_new(analyzer, desired_delay_in_samples);
synthesizer.mixer = gfb_mixer_new(analyzer, synthesizer.delay);

