function scalib = baumgartner2013_calibration(s)
%baumgartner2013_calibration  Calibration of listener-specific sensitivity 
% thresholds to experimental performance
%   Usage: scalib = baumgartner2013_calibration(s)
%
%   Input parameter:
%     s       : strucure containing subject's data. It must include the 
%               fields Obj, pe_exp, and qe_exp, representing the
%               listener's HRTF as SOFA object, the baseline local
%               polar RMS error, and the baseline quadrant error rate,
%               respectively.
%
%   Output parameter:
%     scalib  : strucure containing subject's data with calibrated u
%
%   BAUMGARTNER2013_CALIBRATION returns listener data with
%   listener-specific sensitivity thresholds calibrated by joint
%   optimization of PE and QE to minimize mismatch between experimental
%   and predicted results.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2013_calibration.php

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

kv.latseg = [-20,0,20];

scalib = s;
for ss = 1:length(s)
  
  scalib(ss).u = fminsearch(@(u) evaldist(s(ss),u,kv),s(ss).u,...
    optimset('MaxIter',50,'TolX',0.001)...
    );
  amt_disp([num2str(ss,'%2.0u') ' of ' num2str(length(s),'%2.0u') ' calibrated.'],'progress')

end


end


function [distmetric,qeM,peM] = evaldist(s,u,kv)

if S <= 0
  distmetric = Inf;
  return
end

%% LocaMo
qeM = zeros(length(s),1);
peM = zeros(length(s),1);
for ll = 1:length(s)

  for ii = 1:length(kv.latseg)

    s(ll).sphrtfs{ii} = 0;     % init
    s(ll).p{ii} = 0;        % init

    [s(ll).sphrtfs{ii},polang] = extractsp( kv.latseg(ii),s(ll).Obj );
    [s(ll).p{ii},respangs] = baumgartner2013(...
        s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
        'u',u,'lat',kv.latseg(ii),'polsamp',polang);

    [ qe(ii),pe(ii) ] = baumgartner2013_pmv2ppp( ...
        s(ll).p{ii} , polang , respangs , s(ll).target{ii});

    qeM(ll) = qeM(ll) + qe(ii)*s(ll).Ntargets{ii}/sum([s(ll).Ntargets{:}]);
    peM(ll) = peM(ll) + pe(ii)*s(ll).Ntargets{ii}/sum([s(ll).Ntargets{:}]);

  end

  dQE(ll) = s(ll).qe_exp - qeM(ll);
  dPE(ll) = s(ll).pe_exp - peM(ll);

end

[qe_chance,pe_chance] = baumgartner2013_pmv2ppp(ones(49,44));
distmetric =  (dQE/qe_chance).^2 + (dPE/pe_chance).^2; % Joint distance metric of QE and PE (standardized scatter)

end
