function definput=arg_baumgartner2014_pmv2ppp(definput)

definput.keyvals.p=ones(72,44);
definput.keyvals.rang=-90:5:269;
definput.keyvals.tang=[-30:5:70,80,100,110:5:210];
definput.keyvals.exptang=[];

definput.flags.print = {'noprint','print'};
definput.flags.chance = {'','chance'};
definput.flags.ppp = {'','QE_PE_EB','QE','PE','EB','absPE'};


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_baumgartner2014_pmv2ppp.php



