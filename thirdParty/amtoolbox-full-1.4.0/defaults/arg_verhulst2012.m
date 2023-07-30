function definput=arg_verhulst2012(definput)
% ARG_VERHULST2012
%
%   #License: GPL
%   #Author: Piotr Majdak (2021)
%% General
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_verhulst2012.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
definput.flags.disp = {'no_debug','debug'};

definput.keyvals.normalize = []; % leave empty for no normalization, or provide a binary vector enabling normalization at each channel. 
definput.keyvals.subject = 1; % standard subject controls the cochlear irregulatiries.
definput.keyvals.irr = []; % leave empty for irregularities in all channels or provide a binary vector enabling irregularities at each channel.




