function [prior,mrs,s] = baumgartner2016_parametrization(s)
% baumgartner2016_parametrization - Joint optimization of model parameters
%   Usage: [gamma,prior] = baumgartner2016_parametrization(s,kv)
%
%   Input parameter:
%     s       : strucure containing subject's data. It must include the 
%               following fields:
%               Obj ... the listener's HRTF as SOFA object. 
%               itemlist ... the listener's response patterns. (See help
%               localizationerror)
%
%   Output parameter:
%     gamma   : degree of selectivity in 1/dB
%
%     prior   : prior expectation paramter
%
%   BAUMGARTNER2016_PARAMETRIZATION(...) jointly optimizes the degree of 
%   selectivity Gamma, the response scatter epsilon induced by
%   sensorimotor mapping, and the listener-specific sensitivity S_l.
%
%   Examples:
%   ---------
%
%   This example shows how to parametrize the model according to the data
%   from baumgartner2016 :
%     
%     s = data_baumgartner2016('Long','model'); % Load the experimental data
%     [gamma,prior] = baumgartner2016_parametrization(s);
%
%   See also: baumgartner2016, data_baumgartner2016
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling the effects of
%     sensorineural hearing loss on auditory localization in the median
%     plane. Trends in Hearing, 20:1-11, 2016.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2016_parametrization.php

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

%% Evaluate performance as a function of the lateral angle

definput=arg_baumgartner2016;

disp('Parametrization only done for median plane!')
latseg = 0;%-60:20:60; % centers of lateral segments
dlat =  30;  % lateral range (+-) of each segment

for ll = 1:length(s)

  s(ll).target = [];
  s(ll).response = [];
  s(ll).Nt = [];
  s(ll).baseline.pe_exp = real(localizationerror(s(ll).itemlist,'rmsPmedianlocal'));
  s(ll).baseline.qe_exp = real(localizationerror(s(ll).itemlist,'querrMiddlebrooks'));
  s(ll).pe_exp_lat = zeros(1,length(latseg));
  s(ll).qe_exp_lat = zeros(1,length(latseg));
  for ii = 1:length(latseg)

    latresp = s(ll).itemlist(:,7);
    idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
    s(ll).mm2 = s(ll).itemlist(idlat,:);

    s(ll).mm2(:,7) = 0; % set lateral angle to 0deg such that localizationerror works also outside +-30deg

    s(ll).pe_exp_lat(ii) = real(localizationerror(s(ll).mm2,'rmsPmedianlocal'));
    s(ll).qe_exp_lat(ii) = real(localizationerror(s(ll).mm2,'querrMiddlebrooks'));

    s(ll).target{ii} = real(s(ll).mm2(:,6)); % polar angle of target
    s(ll).response{ii} = real(s(ll).mm2(:,8)); % polar angle of response
    s(ll).Nt{ii} = length(s(ll).target{ii});

  end
end

%% Optimize

x0 = [0,25]; % init of Gamma resp. epsilon
xopt = fminsearch(@(x) evaldist_baumgartner2016parametrization(s,x,latseg,definput),x0,...
    optimset('Display','iter','MaxIter',50,'TolX',1,'PlotFcns',@optimplotx)...
    );
  
% xopt = fminbnd(@(x) evaldist_baumgartner2016parametrization(s,x,latseg,definput),eps,100,...
%   optimset('Display','iter','MaxIter',50,'TolX',0.1,'PlotFcns',@optimplotx)...
%   );

gamma = definput.keyvals.gamma
prior = xopt(1)%definput.keyvals.prior;
mrs = xopt(2)


end


function [distmetric] = evaldist_baumgartner2016parametrization(s,x,latseg,definput)

gamma = definput.keyvals.gamma;
prior = x(1); % x(2);
epsilon = x(2);%17; 

%% Calibrate the sensitivity
kv.mrsmsp = epsilon;
kv.gamma = gamma;
kv.do = 1;
c.latseg = 0;%[-20,0,20];
% c.SPL = [s.SPL];
% for ii = 1:length(s)
%   c.stim{ii} = s(ii).stim;
% end
% s = baumgartner2016_calibration(s,kv,c,0.5);
s = baumgartner2016_calibration(s,kv);

%% Total number of targets
Nt = zeros(length(s),1); %init
for ll = 1:length(s)
  Nt(ll) = sum([s(ll).Nt{:}]);
end
Ntotal = sum(Nt);

%% LocaMo
dQEsq = 0;
dPEsq = 0;
for ll = 1:length(s)

  Nt(ll) = sum([s(ll).Nt{:}]);
  
  for ii = 1:length(latseg)

    if s(ll).Nt{ii} > 0
      s(ll).sphrtfs{ii} = 0;     % init
      s(ll).p{ii} = 0;        % init

      [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
%       [s(ll).p{ii},respangs] = baumgartner2016(...
%           s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).id,s(ll).fs,...
%           'lat',latseg(ii),'polsamp',polang,...
%           'S',s(ll).S,...
%           'mrsmsp',epsilon,'gamma',gamma,'prior',prior,...
%           'SPL',s(ll).SPL,'stim',s(ll).stim,'fsstim',s(ll).fsstim);
      [s(ll).p{ii},respangs] = baumgartner2016(...
          s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).id,s(ll).fs,...
          'lat',latseg(ii),'polsamp',polang,...
          'S',s(ll).S,...
          'mrsmsp',epsilon,'gamma',gamma,'prior',prior);

      [ qe,pe ] = baumgartner2014_pmv2ppp( ...
          s(ll).p{ii} , polang , respangs , s(ll).target{ii});

      dQE_lat = qe - s(ll).qe_exp_lat(ii);
      dPE_lat = pe - s(ll).pe_exp_lat(ii);

      % Accumulate squared errors weighted by number of targets
      dQEsq = dQEsq + dQE_lat.^2 * s(ll).Nt{ii}/Ntotal; 
      dPEsq = dPEsq + dPE_lat.^2 * s(ll).Nt{ii}/Ntotal;
    end

  end

end

% [qe_chance,pe_chance] = pmv2ppp(ones(49,44));
% distmetric =  (dQEsq/qe_chance^2) + (dPEsq/pe_chance^2); % Joint distance metric of QE and PE (normalized by chance performance)

QEmax = 100;
PEmax = 90;
distmetric =  sqrt((dQE/QEmax).^2 + (dPE/PEmax).^2);

end
