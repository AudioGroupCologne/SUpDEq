function [output, analyzer] = gfb_analyzer_process(analyzer, input)
%GFB_ANALYZER_PROCESS  Process input data
%   Usage: [output, analyzer] = gfb_analyzer_process(analyzer, input);
%
%   Input parameters:
%     analyzer : A gfb_analyzer struct created from gfb_analyzer_new. The
%                analyzer will be returned (with updated filter states) as
%                the second output parameter
%     input    : Either a row vector containing the input signal to process, or
%                a matrix containing different input signals for the different
%                bands.  Different rows correspond to different filter bands,
%                while different colums correspond to different times. 
%
%   Output parameters:
%     output   : A matrix containing the analyzer's output signals.  The
%                rows correspond to the filter bands
%
%   The analyzer processes the input data.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/gfb_analyzer_process.php

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

% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Sep 2003, Nov 2006, Jan 2007

warning('Warning: GFB_ANALYZER_PROCESS will be removed in a future release. Use hohmann2002_process instead. ');


if (analyzer.fast)
  % use matlab extension for fast computation.
  [output, analyzer] = gfb_analyzer_fprocess(analyzer, input);
else
  number_of_bands = length(analyzer.center_frequencies_hz);
  output = zeros(number_of_bands, length(input));
  for band = [1:number_of_bands]
    [output(band,:), analyzer.filters(band)] = ...
        gfb_filter_process(analyzer.filters(band), ...
                           input );
  end
end

