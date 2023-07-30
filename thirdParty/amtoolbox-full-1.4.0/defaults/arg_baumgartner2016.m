function definput=arg_baumgartner2016(definput)

definput = arg_baumgartner2016_calibration(definput);
definput = arg_localizationerror(definput);
definput = arg_baumgartner2014_pmv2ppp(definput);
definput = arg_amt_cache(definput);

definput.keyvals.ID = 'NHx';
definput.keyvals.Condition = 'baseline';%'Long';

tmp = amt_load('baumgartner2016','temstim.mat'); % frozen version of noise(8e3,1,'white'); potential alternative: mls(13)
definput.keyvals.temstim = tmp.temstim;

definput.flags.fbank = {'zilany2014','zilany2007','gammatone'};   % ,'lopezpoveda2001'
definput.keyvals.fiberTypes=1:3;  % IHC scaling factor: 1 denotes normal IHC function; 0 denotes complete IHC dysfunction.
definput.flags.redoSpectralAnalysis= {'','redoSpectralAnalysis'};
definput.flags.SPLtemAdapt = {'','SPLtemAdapt'};
definput.flags.normalHearingTemplate = {'','NHtem'};

definput.flags.Ifw = {'nointensityweighting','intensityweighting'};
definput.flags.diff = {'','diff','nodiff'};
definput.flags.sensitivitymapping = {'sigmoid','exp','Gauss','threshold'};
definput.flags.comparisonprocess = {'isd','corr'};
definput.flags.featureextractor = {'psge','dcn'};
definput.flags.fibertypeseparation = {'ftcum','ftopt',''};

definput.flags.headphonefilter = {'','headphone'};
definput.flags.middleearfilter = {'no_middleear','middleear'};
definput.flags.ihc = {'noihc','ihc'};    
definput.flags.regularization = {'regular','noregular'};
definput.flags.motoricresponsescatter = {'mrs','nomrs'};

definput.flags.settings = {'notprint','print'};

definput.keyvals.fs=48000;      % Hz
definput.keyvals.S=0.5;         % listener-specific sensitivity parameter
definput.keyvals.lat=0;         % deg
definput.keyvals.stim=[];
definput.keyvals.fsstim=48e3;
definput.keyvals.space=1;       % No. of ERBs (Cams) 
definput.keyvals.do=1;
definput.keyvals.flow=700;      % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.SPL = 60; 	% dB SPL
definput.keyvals.SPLtem = [40,80];	% dB SPL
definput.keyvals.SL = [];       % db/ERB; spectral density of target sound re absolut detection threshold
definput.keyvals.bwcoef=13;     % steepness in degrees of binaural weighting function
definput.keyvals.polsamp=[-30:5:70 80 100 110:5:210];  % polar sampling (for regular)
definput.keyvals.rangsamp=5;    % equi-polar sampling of response angles
definput.keyvals.mrsmsp=17;%25.5;     % degrees
definput.keyvals.gamma=6;       % slope of psychometric function
definput.keyvals.prior=0;       % parameter of Pratt prior
definput.keyvals.priordist.x = [-90,270]; % angles of prior distribution
definput.keyvals.priordist.y = [1,1];     % values of prior distribution

definput.keyvals.SimDL=eps; % Difference limen of similarity estimation
definput.keyvals.SimThresh=[]; % Threshold of similarity estimation

definput.keyvals.numCF = 28;      % # CFs for zilany model
definput.keyvals.fsmod = 100e3; % Hz, sampling rate of zilany model
definput.keyvals.cohc=1;      % OHC scaling factor: 1 denotes normal OHC function; 0 denotes complete OHC dysfunction.
definput.keyvals.cihc=1;      % IHC scaling factor: 1 denotes normal IHC function; 0 denotes complete IHC dysfunction.

definput.keyvals.gammashortfact = 1; % adaptation factor for duration-dependent gamma
definput.keyvals.Sshortfact = 1; % adaptation factor for duration-dependent sensitivity
definput.keyvals.psgeshort = 1; % inhibitory strength for short sounds

definput.keyvals.GT_minSPL = 20; % minimum representable SPL
definput.keyvals.GT_maxSPL = 110; % maximum representable SPL

definput.keyvals.tiwin = inf; % temporal integration window in seconds

definput.keyvals.mgs = 6; % maximum gradient sensitivity


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_baumgartner2016.php



