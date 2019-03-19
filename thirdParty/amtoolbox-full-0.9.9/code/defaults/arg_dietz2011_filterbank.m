function definput=arg_dietz2011_filterbank(definput)
 
  % Parameters for filtering the haircell output
  definput.keyvals.filter_order = 2;            % used for both env and fine
  definput.keyvals.filter_attenuation_db = 10;  % used for both env and fine
  % Finestructure filter
  definput.keyvals.fine_filter_finesse = 3;
  % Envelope/modulation filter
  definput.keyvals.mod_center_frequency_hz = 135;
  definput.keyvals.mod_filter_finesse = 8; % => bandwidth: 16.9 Hz
  % ILD filter
  definput.keyvals.level_filter_order = 2;
  definput.keyvals.level_filter_cutoff_hz = 30;

%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/defaults/arg_dietz2011_filterbank.php

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

