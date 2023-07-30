function definput = arg_hauth2020(definput)
% ARG_HAUTH2020
%
%   #License: GPL
%   #Author: Piotr Majdak (2021)
%% General
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_hauth2020.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
definput.flags.disp = {'no_debug','debug'};
definput.keyvals.long_term = 1; % use long-term model (1) or short-term model (0)
definput.flags.timescale = {'longterm', 'shortterm'};

%% Filterbank
definput.keyvals.Filterbank = 'GT';% define filtering: GT is gammatone filter
definput.keyvals.fmin = 150;
definput.keyvals.fmax = 8500; % highest filter of gammatone filterbank
definput.keyvals.f_target = 500; % specified frequeny that will have a matched filter
%definput.keyvals.bin_err = 1; % enable(1)/disable(0) binaural processing inaccuracies (important for Monte-Carlo Simulations)
definput.keyvals.ERB_factor = 1;% define bandwidth of filters
definput.keyvals.OptSigs = [];
definput.flags.bin_err = {'binauralinaccuracies', 'no_binauralinaccuracies'};


