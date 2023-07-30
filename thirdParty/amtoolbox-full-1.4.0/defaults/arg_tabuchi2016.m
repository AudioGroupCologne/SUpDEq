function definput=arg_tabuchi2016(definput)
% ARG_TABUCHI2016
%
%   #License: GPL
%   #Author: Hisaaki Tabuchi (2022)
%   #Author: Clara Hollomey (2023)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_tabuchi2016.php


definput.keyvals.ListenerStr = 'GrandMeanSd'; % to get the saved data
definput.keyvals.n1 = 34;
definput.keyvals.n2 = 46;
definput.keyvals.Cvec = -1.5:0.25:1; % This range of Cvec was tested by considering the presumed curvature, -0.5. See the footnote 2 in Tabuchi et al. (2016). 
definput.keyvals.GmaxVec = 0:70;
definput.keyvals.x_vec = 0:0.1:110; % target level, dB SPL with 0.1 dB step
definput.keyvals.CondStr_vec = {'OffFreqPre', 'OnFreqPre'};
definput.keyvals.mlvl_vec = [60 90];
definput.keyvals.GmaxToGetK = 34;
definput.keyvals.gammaVec = [60]; % for IHC 
definput.keyvals.beta = 1; % for IHC 
