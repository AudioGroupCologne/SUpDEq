function fb = hohmann2002(fs,flow,basef,fhigh,filters_per_ERBaud,varargin)
%HOHMANN2002  Invertible Gammatone filterbank
%   Usage:  fb = hohman2002(fs, flow, basef, fhigh, filters_per_ERBaud)
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
%                         constant bandwidth factor. Default is 1.0.
%
%     'l'                 Scaling factor l from Eq. 13 in Hohmann (2002).
%                         Default is 24.7.
%  
%     'q'                 Scaling factor q from Eq. 13 in Hohmann (2002).
%                         Default is 9.265.
%  
%   References:
%     V. Hohmann. Frequency analysis and synthesis using a gammatone
%     filterbank. Acta Acustica united with Acoustica, 88(3):433--442, 2002.
%     
%
%   See also: demo_hohmann2002, demo_gammatone, hohmann2002_mixer,
%             dietz2011_filterbank, hohmann2002_synth,
%             hohmann2002_filter, hohmann2002_freqz,
%             hohmann2002_delay, hohmann2002_clearstate, hohmann2002_process, 
%             kelvasa2015_localize, exp_hohmann2002, exp_gammatone, hauth2020
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/hohmann2002.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified  
%   #Author   : Universitaet Oldenburg, tp (2002 - 2007)
%   #Author   : Piotr Majdak (2016)
%   Adapted from function gfb_analyzer_new

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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



