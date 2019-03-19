function definput=arg_baumgartner2014(definput)
  
definput.flags.regularization = {'regular','noregular'};
definput.flags.motoricresponsescatter = {'mrs','nomrs'};
definput.flags.settings = {'notprint','print'};

% CP-Falgs:
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/defaults/arg_baumgartner2014.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
definput.flags.cp={'fwstd','std','xcorr'};

definput.keyvals.fs=48000;      % sampling rate in Hz
definput.keyvals.S=0.5;         % listener-specific sensitivity
definput.keyvals.lat=0;         % lateral angle in deg
definput.keyvals.stim=[];
definput.keyvals.fsstim=[];
definput.keyvals.space=1;       % No. of ERBs (Cams) 
definput.keyvals.flow=700;      % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.do=1;          % differential order
definput.keyvals.bwcoef=13;     % binaural weighting coefficient in deg
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

