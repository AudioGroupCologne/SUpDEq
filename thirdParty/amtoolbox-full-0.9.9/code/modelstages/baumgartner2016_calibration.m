function scalib = baumgartner2016_calibration(s,varargin)
%baumgartner2016_calibration  Calibration of listener-specific sensitivity thresholds to experimental performance
%   Usage: scalib = baumgartner2016_calibration(s)
%
%   Input parameter:
%     s       : strucure containing subject's data. It must include the 
%               fields Obj, baseline.pe_exp, and baseline.qe_exp, representing 
%               the listener's HRTF as SOFA object, the baseline local
%               polar RMS error, and the baseline quadrant error rate,
%               respectively.
%
%   Output parameter:
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
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2016_calibration.php

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

% AUTHOR : Robert Baumgartner

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
    xopt = fminsearch(@(x) evaldist(s(ss),x,kv,flags,c.SPL(ss),c.stim{ss}),...
        [mean(kv.Srange) mean(kv.prange)],...
        optimset('MaxIter',kv.MaxIter,'TolX',kv.TolX,'Display','iter','TolFun',1e-10)...
        );
    scalib(ss).S = real(xopt(1));
    scalib(ss).prior = real(xopt(2))*100; % to equalize sensitivity between S and prior
  else 
  switch flags.optimization
    
    case 'fminsearch'
      xopt = fminsearch(@(x) evaldist(s(ss),x,kv,flags,c.SPL(ss),c.stim{ss}),mean(kv.Srange),...
        optimset('MaxIter',kv.MaxIter,'TolX',kv.TolX,'Display','iter','TolFun',1e-10)...
        );
      scalib(ss).S = real(xopt);
      
    case 'fminbnd'
      xopt = fminbnd(@(x) evaldist(s(ss),x,kv,flags,c.SPL(ss),c.stim{ss}),kv.Srange(1),kv.Srange(2),...
        optimset('MaxIter',kv.MaxIter,'TolX',kv.TolX,'Display','iter','TolFun',1e-10)...
        );
      scalib(ss).S = real(xopt);
  
    case 'search'
      iS = kv.Srange(1):(kv.Srange(2)-kv.Srange(1))/100:kv.Srange(2);
      distmetric = zeros(length(iS),1);
      for ii = 1:length(iS)
        distmetric(ii) = evaldist(s(ss),iS(ii),kv,flags,c.SPL(ss),c.stim{ss});
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


function [distmetric,qeM,peM] = evaldist(s,x,kv,flags,SPL,stim)

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
