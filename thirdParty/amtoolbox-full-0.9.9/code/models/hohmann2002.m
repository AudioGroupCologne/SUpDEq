function fb = hohmann2002(fs,flow,basef,fhigh,filters_per_ERBaud,varargin)
%HOHMANN2002  Construct a new filterbank within the HOHMANN2002 framework
%   Usage:  fb = hohman2002(fs,flow, basef, fhigh,filters_per_ERBaud)
%
%   Input parameters:
%      fs                 : The sampling frequency of the signals on which
%                           the filterbank will operate
%      flow               : The lowest center frequency of a filter
%      basef              : "base frequency". One of the filters will
%                           have exactly this center frequency. Must be >= flow
%      fhigh              : The highest center frequency of a filter.
%                           Must be >=  basef
%      filters_per_ERBaud : The density of gammatone filters on the ERB
%                           scale.
%
%   Output parameters:
%      fb           : The constructed filterbank object.
%    
%   hohman2002 constructs a new fb object implementing
%   the analysis part of a gammatone filterbank as described
%   in Hohmann (2002).
%
%   fb contains several all-pole gammatone filters; each
%   one with a bandwidth of 1 ERB (optional: times bandwidth_factor),
%   and an order of gamma_order.
%
%   The center frequencies of the individual filters are computed as
%   described in section 3 of Hohmann (2002).
%
%   HOHMANN2002 takes the following flags at the end of the line of input
%   arguments:
%
%     'gamma_order'       The order of the gammatone filters in this
%                         filterbank. Default is 4.
%
%     'bandwidth_factor'  The bandwidth parameter of the individual filters
%                         is calculated from the Equivalent Rectangular
%                         Bandwidth (ERB) according to equation 14 in
%                         Hohmann (2002). ERB is taken from the Glasberg &
%                         Moore formula for a specific center frequency
%                         (equation 13 in Hohmann (2002)).
%                         Using this parameter, it is possible to widen or
%                         narrow all filters of the filterbank with a
%                         constant bandwidth factor. Default is 1.0
%
%     'l'                 Scaling factor l from Eq. 13 in Hohmann (2002).
%                         Default is 24.7.
%  
%     'q'                 Scaling factor q from Eq. 13 in Hohmann (2002).
%                         Default is 9.265.
%  
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/hohmann2002.php

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
  
% author   : Universitaet Oldenburg, tp (Jan 2002, Jan, Sep 2003, Nov 2006, Jan 2007)
% Adapted to AMT (PM, Jan 2016) from function gfb_analyzer_new

definput.keyvals.l=24.7; % see equation (17) in [Hohmann 2002]
definput.keyvals.q=9.265; % see equation (17) in [Hohmann 2002]
definput.keyvals.gamma_order=4;
definput.keyvals.gaincalc_iterations=100; % number of iterations for the approximation of the gain factors
definput.keyvals.bandwidth_factor=1.0;

[~,kv]  = ltfatarghelper({'l','q','gamma_order','gaincalc_iterations','bandwidth_factor'},definput,varargin);

fb.type = 'gfb_analyzer';
fb.fs = fs;
fb.flow = flow;
fb.basef = basef;
fb.fhigh = fhigh;
fb.filters_per_ERBaud = filters_per_ERBaud;
fb.bandwidth_factor = kv.bandwidth_factor;
fb.fast = 0;

fb.L = kv.l;
fb.Q = kv.q;
fb.gamma_order = kv.gamma_order;
fb.gaincalc_iterations = kv.gaincalc_iterations;

fb.center_frequencies_hz = ...
    erbspacebw(flow,fhigh,1/filters_per_ERBaud,basef);

% This loop actually creates the gammatone filters:
for band = 1:length(fb.center_frequencies_hz)
  center_frequency_hz = fb.center_frequencies_hz(band);

  % Construct gammatone filter:
  fb.filters(1,band) = ...
      hohmann2002_filter(fs, center_frequency_hz, kv.gamma_order, kv.bandwidth_factor);
end


