function [output, mixer] = gfb_mixer_process(mixer, input)
% [output, mixer] = gfb_mixer_process(mixer, input)
%
% gfb_mixer_process computes a weighted sum of the different bands present
% in input.
%
%  Input parameters:
% mixer   A gfb_mixer structure as returned by gfb_mixer_new.  The mixer
%         contains the gain factors for the weighted sum.  A copy of mixer
%         will be returned in the second return parameter
% input   an NxM matrix, where N equals the number of gain factors (bands)
%         of the mixer
% output  an 1xM vector containing the weighted sums of each comlumn
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/gfb_mixer_process.php

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

% filename : gfb_mixer_process.m

warning('Warning: GFB_MIXER_PROCESS will be removed in a future release. Use hohmann2002_process instead. ');

output = mixer.gains * input;


%OLDFORMAT

