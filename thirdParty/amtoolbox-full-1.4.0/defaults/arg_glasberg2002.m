function definput = arg_glasberg2002(definput)

definput.keyvals.fs = 32000;
definput.keyvals.flow = 20;
definput.keyvals.fhigh = 16000;

% filter order as in glasberg2002
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_glasberg2002.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
definput.keyvals.order = 4096;
definput.keyvals.fftLen = 2048; % according to glasberg2002
definput.keyvals.vLimitingIndices = [ 4050, 2540, 1250, 500, 80];%f-boundary vector
% compute windows
definput.keyvals.hannLenMs = [2,4,8,16,32,64]; % hanning window size (glasberg2002) in ms
definput.keyvals.timeStep = 0.001; % 1ms steps as in glasberg2002
definput.flags.compensationtype = {'tfOuterMiddle1997','tfOuterMiddle2007','specLoud'};

