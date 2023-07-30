function definput=arg_baumgartner2014(definput)
  
% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_baumgartner2014.php


definput.flags.regularization = {'regular','noregular'};
definput.flags.motoricresponsescatter = {'mrs','nomrs'};
definput.flags.settings = {'notprint','print'};

% CP-Falgs:
definput.flags.cp={'fwstd','std','xcorr'};

definput.keyvals.fs=48000;      % sampling rate in Hz
definput.keyvals.S=0.7;         % listener-specific sensitivity (mean S of data_baumgartner2014 pool)
definput.keyvals.lat=0;         % lateral angle in deg
definput.keyvals.stim=[];
definput.keyvals.fsstim=[];
definput.keyvals.space=1;       % No. of ERBs (Cams) 
definput.keyvals.flow=700;      % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.do=1;          % differential order
definput.keyvals.bwcoef=13;     % binaural weighting coefficient in deg
definput.keyvals.tang=[-30:5:70 80 100 110:5:210];  % polar sampling 
definput.keyvals.polsamp=[-30:5:70 80 100 110:5:210];  % polar sampling (for regular flag)
definput.keyvals.rangsamp=5;    % equi-polar sampling of response angles
definput.keyvals.mrsmsp=17;     % response scatter (in deg) in the median plane induced by non-auditory processes
definput.keyvals.gamma=6;       % degree of selectivity in 1/dB
  

definput.flags.localizationerror = {'','accL','precL','precLcentral','accP','precP','querr','accE',...
  'absaccE','absaccL','accabsL','accPnoquerr','precE','querr90','accE','precPmedian',...
  'precPmedianlocal','precPnoquerr','rmsL','rmsPmedianlocal','rmsPmedian',...
  'querrMiddlebrooks','corrcoefL','corrcoefP','SCC',...
  'gainLstats','gainL','pVeridicalL','precLregress'...
  'sirpMacpherson2000','gainPfront','gainPrear','gainP','perMacpherson2003',...
  'pVeridicalPfront','pVeridicalPrear','precPregressFront','precPregressRear',...
  'QE_PE_EB','QE','PE','EB','absPE'};


