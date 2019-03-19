function [output, synthesizer] = gfb_synthesizer_process(synthesizer, input)
% [output, synthesizer] = gfb_synthesizer_process(synthesizer, input)
%
% The synthesizer will resynthesize the given input.
%
%  Input parameters:
% synthesizer  A synthesizer structure as created by gfb_synthesizer_new. A
%              copy of the synthesizer object with an updated internal state
%              is returned in the second return parameter
% input        A matrix containing the (possibly processed) complex output of
%              the analyzer corresponding to this synthesizer.  The number of
%              rows in input must match the number of filter bands
% output       The synthesized output signal
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/gfb_synthesizer_process.php

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

% filename : gfb_synthesizer_process.m

warning('Warning: GFB_SYNTHESIZER_PROCESS will be removed in a future release. Use hohmann2002_process instead. ');

[output, synthesizer.delay] = gfb_delay_process(synthesizer.delay, input);
[output, synthesizer.mixer] = gfb_mixer_process(synthesizer.mixer, output);

%OLDFORMAT

