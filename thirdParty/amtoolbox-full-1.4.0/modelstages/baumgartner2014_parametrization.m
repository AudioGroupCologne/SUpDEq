function [gamma,epsilon,s] = baumgartner2014_parametrization(s)
% BAUMGARTNER2014_PARAMETRIZATION Joint optimization of model parameters
%   Usage: [gamma,epsilon] = baumgartner2014_parametrization(s)
%
%   Input parameters:
%     s       : strucure containing subject's data. It must include the 
%               following fields:
%               Obj ... the listener's HRTF as SOFA object. 
%               itemlist ... the listener's response patterns. (See help
%               localizationerror)
%
%   Output parameters:
%     gamma   : degree of selectivity in 1/dB
%
%     epsilon : response scatter in degrees induced by sensorimotor mapping
%
%   BAUMGARTNER2014_PARAMETRIZATION(...) jointly optimizes the degree of 
%   selectivity Gamma, the response scatter epsilon induced by
%   sensorimotor mapping, and the listener-specific sensitivity S_l.
%
%   Examples:
%   ---------
%
%   This example shows how to parametrize the model according to the data
%   from baumgartner2014 :
%     
%     s = data_baumgartner2014('baseline'); % Load the experimental data
%     [gamma,epsilon] = baumgartner2014_parametrization(s);
%
%   See also: baumgartner2014, data_baumgartner2014
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2014_parametrization.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA CircStat M-SIGNAL M-Stats O-Statistics
%   #Author: Robert Baumgartner (2014), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Evaluate performance as a function of the lateral angle

latseg = -60:20:60; % centers of lateral segments
dlat =  10;  % lateral range (+-) of each segment

for ll = 1:length(s)

  s(ll).target = [];
  s(ll).response = [];
  s(ll).Nt = [];
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

x0 = [6,17]; % init of Gamma resp. epsilon
xopt = fminsearch(@(x) local_evaldistbaumgartner2014parametrization(s,x),x0,...
    optimset('Display','iter','MaxIter',50,'TolX',1)...
    );

gamma = xopt(1);
epsilon = xopt(2);

end


function [distmetric] = local_evaldistbaumgartner2014parametrization(s,x)

gamma = x(1);
epsilon = x(2);
latseg = -60:20:60; % centers of lateral segments

%% Calibrate the sensitivity
kv.mrsmsp = epsilon;
kv.gamma = gamma;
kv.do = 1;
s = baumgartner2014_calibration(s,kv);

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
      [s(ll).p{ii},respangs] = baumgartner2014(...
          s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).Obj.Data.SamplingRate,...
          'S',s(ll).S,'lat',latseg(ii),'polsamp',polang,...
          'mrsmsp',epsilon,'gamma',gamma);

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

[qe_chance,pe_chance] = baumgartner2014_pmv2ppp(ones(49,44));
distmetric =  (dQEsq/qe_chance^2) + (dPEsq/pe_chance^2); % Joint distance metric of QE and PE (normalized by chance performance)

end


