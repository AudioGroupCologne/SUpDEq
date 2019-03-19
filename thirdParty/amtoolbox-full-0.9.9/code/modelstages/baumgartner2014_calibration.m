function scalib = baumgartner2014_calibration(s,kv,TolX)
%baumgartner2014_calibration - Calibration of the model
%   Usage: scalib = baumgartner2014_calibration(s)
%
%   Input parameter:
%     s       : strucure containing subject's data. It must include the 
%               fields Obj, pe_exp, and qe_exp, representing the
%               listener's HRTF as SOFA object, the baseline local
%               polar RMS error, and the baseline quadrant error rate,
%               respectively. Optionally, the structure can include the
%               field target, a cell array representing the polar target
%               angles for each lateral segment.
%
%   Output parameter:
%     scalib  : strucure containing subject's data with calibrated S
%
%   BAUMGARTNER2014_CALIBRATION returns listener data with
%   listener-specific sensitivity thresholds calibrated by joint
%   optimization of PE and QE to minimize mismatch between experimental
%   and predicted results.
%
%   BAUMGARTNER2014_CALIBRATION accepts the following optional parameters:
%
%     kv        Key-value pairs according to baumgartner2014
%
%     TolX      Optimization tolerance of listener-specific sensitivity.
%               Default is 1e-3.
%
%   See also: baumgartner2014
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2014_calibration.php

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

if not(exist('kv','var'))
  definput.import={'baumgartner2014'};
  [~,kv]=ltfatarghelper({},definput,{});
end

if not(isfield(kv,'latseg'))
  kv.latseg = [-20,0,20];
end

if not(isfield(s,'target'))
  amt_disp('Calibration accuracy could be enhanced by providing the target polar-angles.')
end

if not(exist('TolX','var'))
  TolX = 0.001;
end
 
scalib = s;
for ss = 1:length(s)
  
  scalib(ss).S = fminsearch(@(S) evaldist(s(ss),S,kv),s(ss).S,...
    optimset('MaxIter',50,'TolX',TolX)...
    );
%   [~,scalib(ss).qe_pred,scalib(ss).pe_pred] = evaldist(s(ss),S,kv);
  amt_disp([num2str(ss,'%2.0u') ' of ' num2str(length(s),'%2.0u') ' calibrated.'],'progress')

end


end


function [distmetric,qeM,peM] = evaldist(s,S,kv)

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
    [s(ll).p{ii},respangs] = baumgartner2014(...
        s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).Obj.Data.SamplingRate,...
        'S',S,'lat',kv.latseg(ii),'polsamp',polang,...
        'mrsmsp',kv.mrsmsp,'gamma',kv.gamma,'do',kv.do);

    if isfield(s,'target')
      
      [ qe(ii),pe(ii) ] = baumgartner2014_pmv2ppp( ...
          s(ll).p{ii} , polang , respangs , s(ll).target{ii});
      latweight = length(s(ll).target{ii})/length(cat(1,s(ll).target{:}));
      qeM(ll) = qeM(ll) + qe(ii)*latweight;
      peM(ll) = peM(ll) + pe(ii)*latweight;
      
    else
      
      [ qe(ii),pe(ii) ] = baumgartner2014_pmv2ppp( ...
          s(ll).p{ii} , polang , respangs);
      qeM(ll) = qeM(ll) + qe(ii)*1/length(kv.latseg);
      peM(ll) = peM(ll) + pe(ii)*1/length(kv.latseg);
      
    end

  end

  dQE(ll) = s(ll).qe_exp - qeM(ll);
  dPE(ll) = s(ll).pe_exp - peM(ll);

end

[qe_chance,pe_chance] = baumgartner2014_pmv2ppp(ones(49,44));
distmetric =  (dQE/qe_chance).^2 + (dPE/pe_chance).^2; % Joint distance metric of QE and PE (standardized scatter)

end
