function scalib = baumgartner2016_calibration(s,varargin)
%BAUMGARTNER2016_CALIBRATION  Calibration of listener-specific sensitivity thresholds to experimental performance
%   Usage: scalib = baumgartner2016_calibration(s)
%
%   Input parameters:
%     s       : strucure containing subject's data. It must include the 
%               fields Obj, baseline.pe_exp, and baseline.qe_exp, representing 
%               the listener's HRTF as SOFA object, the baseline local
%               polar RMS error, and the baseline quadrant error rate,
%               respectively.
%
%   Output parameters:
%     scalib  : strucure containing subject's data with calibrated u
%
%   BAUMGARTNER2016_CALIBRATION returns listener data with
%   listener-specific sensitivity thresholds calibrated by joint
%   optimization of PE and QE to minimize mismatch between experimental
%   and predicted results.
%
%   BAUMGARTNER2016_CALIBRATION accepts the following optional parameters:
%
%     'Srange',Sr     Define the sensitivity range. Default is [eps,10].
%
%     'prange',pr     Define the prior range. Default is [0,1].
%
%     'latseg',ls     Define lateral segment(s) of data used for
%                     calibration. Default value is 0 deg.
%
%     'c',c           Structure for optional definition of listener-specific 
%                     settings like playback level or stimulus.
%
%     'TolX',tx       Minimum tolerance of optimization argument (see help
%                     optimset for details).
%
%     'MaxIter',mi    Maximum number of optimization iterations (see help
%                     optimset for details).
%
%   BAUMGARTNER2016_CALIBRATION accepts the following flags:
%
%     'calibprior'    Calibrate also expectation prior.
%
%     'fminbnd'       Use fminbnd routine for calibration. This is the
%                     default.
%
%     'fminsearch'    Use fminsearch routine for calibration.
%
%     'search'        Try all possibilities.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2016_calibration.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA M-Signal M-Stats O-Statistics
%   #Author: Robert Baumgartner (2016), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.import={'baumgartner2016'};
definput.keyvals.c = {};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);

c = kv.c;
if isempty(c)
  c.SPL = kv.SPL*ones(length(s),1);
  c.stim = cell(length(s),1);
end

if flags.do_gammatone && kv.tiwin > .5
  kv.Srange = [-1,5];
elseif flags.do_gammatone && kv.tiwin < .5
  kv.Srange = [-5,10];
end

scalib = s;
for ss = 1:length(s)
  if flags.do_calibprior
    xopt = fminsearch(@(x) local_evaldist(s(ss),x,kv,flags,c.SPL(ss),c.stim{ss}),...
        [mean(kv.Srange) mean(kv.prange)],...
        optimset('MaxIter',kv.MaxIter,'TolX',kv.TolX,'Display','iter','TolFun',1e-10)...
        );
    scalib(ss).S = real(xopt(1));
    scalib(ss).prior = real(xopt(2))*100; % to equalize sensitivity between S and prior
  else 
  switch flags.optimization
    
    case 'fminsearch'
      xopt = fminsearch(@(x) local_evaldist(s(ss),x,kv,flags,c.SPL(ss),c.stim{ss}),mean(kv.Srange),...
        optimset('MaxIter',kv.MaxIter,'TolX',kv.TolX,'Display','iter','TolFun',1e-10)...
        );
      scalib(ss).S = real(xopt);
      
    case 'fminbnd'
      xopt = fminbnd(@(x) local_evaldist(s(ss),x,kv,flags,c.SPL(ss),c.stim{ss}),kv.Srange(1),kv.Srange(2),...
        optimset('MaxIter',kv.MaxIter,'TolX',kv.TolX,'Display','iter','TolFun',1e-10)...
        );
      scalib(ss).S = real(xopt);
  
    case 'search'
      iS = kv.Srange(1):(kv.Srange(2)-kv.Srange(1))/100:kv.Srange(2);
      distmetric = zeros(length(iS),1);
      for ii = 1:length(iS)
        distmetric(ii) = local_evaldist(s(ss),iS(ii),kv,flags,c.SPL(ss),c.stim{ss});
      end
      [~,Imin] = min(distmetric);
      scalib(ss).S = iS(Imin);
      
      figure; plot(iS,distmetric); title(s(ss).id)
      pause(0.5)
  end
  end
  disp([num2str(ss,'%2.0u') ' of ' num2str(length(s),'%2.0u') ' calibrated.'])

end


end


function [distmetric,qeM,peM] = local_evaldist(s,x,kv,flags,SPL,stim)

% if x(1) < 0 || x(2) <= 0
%   distmetric = Inf;
%   return
% end

% if x < 0
%   distmetric = Inf;
%   return
% end

% SimDL = x(1);
S = real(x(1));
if length(x) == 2
  prior = x(2);
else
  prior = 0;
end

if prior < 0
  distmetric = Inf;
  return
end

%% LocaMo
qeM = zeros(length(s),1);
peM = zeros(length(s),1);
for ll = 1:length(s)
  
  if not(isfield(s,'fsstim'))
    s(ll).fsstim = s(ll).fs;
  end

  for ii = 1:length(kv.latseg)

    s(ll).p{ii} = 0;        % init

    [s(ll).p{ii},respangs,polang] = baumgartner2016(...
        s(ll).Obj,s(ll).Obj,'argimport',flags,kv,...
        'ID',s(ll).id,'fs',s(ll).fs,'mrsmsp',s(ll).mrs,'S',S,'SPL',SPL,...
        'stim',stim,'fsstim',s(ll).fsstim,'priordist',s(ll).priordist);

    [ qe(ii),pe(ii) ] = baumgartner2014_pmv2ppp( ...
        s(ll).p{ii} , polang , respangs , s(ll).target{ii});

    latweight = length(s(ll).target{ii})/length(cat(1,s(ll).target{:}));
    qeM(ll) = qeM(ll) + qe(ii)*latweight;
    peM(ll) = peM(ll) + pe(ii)*latweight;

  end

  dQE(ll) = s(ll).qe_exp - qeM(ll); % s(ll).baseline.qe_exp 
  dPE(ll) = s(ll).pe_exp - peM(ll); % s(ll).baseline.pe_exp 

end

% [qe_chance,pe_chance] = baumgartner2014_pmv2ppp('chance');
% distmetric =  (dQE/qe_chance).^2 + (dPE/pe_chance).^2; % Joint distance metric of QE and PE (standardized scatter)

QEmax = 100;
PEmax = 90;
distmetric =  sqrt((dQE/QEmax).^2 + (dPE/PEmax).^2);

end


