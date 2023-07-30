function definput = arg_chen2011(definput)

definput.flags.outerearcorrection = {'FreeField', 'PDR10', 'Eardrum'};%OuterEarOpt
definput.keyvals.HLohcdB0 = 0;%no hearing impairment
definput.keyvals.HLihcdB0 = 0;
definput.keyvals.HLcf = 0;
%these values are for calculation of the reference cf...these are the CF of
%the auditory filters, not (necessarily) those in inputF
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_chen2011.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
definput.keyvals.cambin = 0.25;
definput.keyvals.flow = 40;
definput.keyvals.fhigh = 17000;


