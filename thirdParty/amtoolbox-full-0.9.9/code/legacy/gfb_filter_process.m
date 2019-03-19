function [output, filter_obj] = gfb_filter_process(filter_obj, input)
%GFB_FILTER_PROCESS  Filter input data
%   Usage: [output, filter] = gfb_filter_process(filter, input)
%
%
%   Input parameters:
%     filter  : A gfb_Filter struct created from gfb_filter_new.  The filter
%               will be returned with an updated filter state as the second
%               return parameter
%     input   : A vector containing the input signal to process
%
%   Output parameters:
%     output  : A vector containing the filter's output signal
%
%   The filter processes the input data.
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/gfb_filter_process.php

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
% date     : Jan 2002, Nov 2006, Jan 2007

warning('Warning: GFB_FILTER_PROCESS will be removed in a future release. Use hohmann2002_process instead. ');

factor = filter_obj.normalization_factor;

% for compatibility of the filter state with the MEX extension, we
% have to multiply the filter state with the filter coefficient before the
% call to filter:
filter_state = filter_obj.state * filter_obj.coefficient;

for i = [1:filter_obj.gamma_order]
  [input, filter_state(i)] = ...
      filter(factor, [1, -filter_obj.coefficient], ...
             input, filter_state(i));
  factor = 1;
end

output = input;

% for compatibility of the filter state with the MEX extension, we
% have to divide the filter state by the filter coefficient after the
% call to filter:
filter_obj.state = filter_state / filter_obj.coefficient;

