function definput = arg_baumgartner2021(definput)

definput.keyvals.tempWin = 1;%0.02; % temporal integration window in sec
definput.keyvals.reflectionOnsetTime = [];
definput.flags.normalize = {'regular','normalize'};
definput.flags.cueProcessing = {'misc','intraaural','interaural'};
definput.flags.lateralInconsistency = {'noLateralInconsistency','lateralInconsistency'};
definput.flags.middleEarFilter = {'','middleEarFilter'};
definput.flags.spectralCueEchoSuppression = {'','spectralCueEchoSuppression'};
definput.keyvals.ILD_JND = 1; % ILD JND from reference
definput.keyvals.ITD_JND = 20e-6; % ITD JND from reference
definput.keyvals.range = 1; % scaling of externalization score
definput.keyvals.offset = 0; % offset of externalization score
definput.flags.decisionStatistics = {'','dprime'};
definput.flags.gradients={'positive','negative','both'};
definput.flags.strategy={'fixedWeights','LTA','MTA','WTA','WSM'};
definput.keyvals.cueWeights = [0.6,0.4];
definput.keyvals.S = [0.20,0.51]; % results taken from exp_baumgartner2020('fig3') 


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_baumgartner2021.php



