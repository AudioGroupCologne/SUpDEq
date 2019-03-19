function obj = hohmann2002_clearstate(obj)
%hohmann2002_clearstate  Reset states of hohmann2002 filters
%   Usage: filter = hohmann2002_clearstate(filter)
%          fb = hohmann2002_clearstate(fb)
%
%   filter = HOHMANN2002_CLEARSTATE(filter) resets the states of the
%   filter created by HOHMANN2002_FILTER
%
%   fb = HOHMANN2002_CLEARSTATE(fb) resets the states of the
%   filterbank fb created by HOHMANN2002
% 
%   delay = HOHMANN2002_CLEARSTATE(delay) resets the states of the
%   delay created by HOHMANN2002_DELAY
%
%   synth = HOHMANN2002_CLEARSTATE(synth) resets the states of the
%   sinthesis filterbank created by HOHMANN2002_SYNTH
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/hohmann2002_clearstate.php

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
  
% Original author: Universitaet Oldenburg, tp (Jan 2002, Nov 2006, Feb 2007)
% Adapted to AMT (PM, Jan 2016) from functions gfb_*_clear_state

if ~isfield(obj,'type'), error('Type of the object missing'); end
switch(obj.type)
  case 'gfb_Filter'
    obj.state = zeros(1, obj.gamma_order);
  case 'gfb_analyzer'
    for band = 1:length(obj.center_frequencies_hz)
        obj.filters(1, band).state = zeros(1, obj.filters(1, band).gamma_order);
    end
  case 'gfb_Delay'
    obj.memory = zeros(size(obj.memory));
  case 'gfb_Synthesizer'
    obj.delay.memory = zeros(size(obj.delay.memory));
  otherwise
    error('Unknown type of HOHMANN2002 filter object');    
end

