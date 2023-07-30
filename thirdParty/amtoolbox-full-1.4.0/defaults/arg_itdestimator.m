function definput=arg_itdestimator(definput)
  
definput.flags.disp = {'no_debug','debug'};
definput.flags.mode = {'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe',...
    'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD'};
definput.flags.lp = {'lp','bb'};
definput.flags.peak = {'hp','fp'};
definput.flags.toaguess = {'noguess','guesstoa'};

definput.keyvals.threshlvl = -10;
definput.keyvals.butterpoly = 10;
definput.keyvals.upper_cutfreq = 3000;
definput.keyvals.lower_cutfreq = 1000;
definput.keyvals.avgtoa = 45;
definput.keyvals.fs = [];

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_itdestimator.php

definput.keyvals.fs = [];

