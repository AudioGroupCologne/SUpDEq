function definput = arg_baumgartner2016_calibration( definput )

definput.keyvals.TolX = 0.005;
definput.keyvals.MaxIter = 100;
definput.keyvals.Srange = [-1,2.5];%[eps,10];
definput.keyvals.prange = [0,1];
definput.keyvals.latseg = 0;
definput.keyvals.dlat = 30;
definput.keyvals.c = {};

definput.flags.recalib={'','recalib'};
definput.flags.prior = {'','calibprior'};
definput.flags.optimization = {'fminbnd','fminsearch','search'};


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_baumgartner2016_calibration.php



