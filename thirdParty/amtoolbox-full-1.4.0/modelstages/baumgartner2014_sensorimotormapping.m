function [ri,rang] = baumgartner2014_sensorimotormapping(si,varargin)
%BAUMGARTNER2014_SENSORIMOTORMAPPING Response scatter induced by localization task
%   Usage:     [ri,rang] = baumgartner2014_sensorimotormapping(si)
%
%   Input parameters:
%     ri      : response index
%
%   Output parameters:
%     si      : similarity index
%     rang    : response polar angles
%
%   BAUMGARTNER2014_SENSORIMOTORMAPPING(...) performs polar-angle
%   interpolation and emulates task-induced response scatter.
%
%   BAUMGARTNER2014_SENSORIMOTORMAPPING accepts the following optional parameters:
%
%     'polsamp',ps   Define the polar-angle sampling of the acoustic data
%                    provided for the current sagittal plane. As default the 
%                    sampling of ARI's HRTFs in the median SP is used, i.e.,
%                    ps = [-30:5:70,80,100,110:5:210] degrees.
%
%     'rangsamp',rs  Define the equi-polar sampling of the response predictions.
%                    The default is rs = 5 degrees.
%
%     'mrsmsp',eps   Set the motoric response scatter eps within the median 
%                    sagittal plane. Default value is 17 degrees.
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2014_sensorimotormapping.php


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


definput.import={'baumgartner2014'};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% Interpolation (regularize polar angular sampling)
if flags.do_regular
    rang0 = ceil(min(kv.polsamp)*1/kv.rangsamp)*kv.rangsamp;    % ceil to kv.rangsamp deg
    rang = rang0:kv.rangsamp:max(kv.polsamp);
    siint = zeros(length(rang),size(si,2));
    for tt = 1:size(si,2)
        siint(:,tt) = interp1(kv.polsamp,si(:,tt),rang,'spline');
    end
    si = siint;
    si(si<0) = 0; % SIs must be positive (necessary due to spline interp)
else
    rang = kv.polsamp;
end


%% Sensorimotor mapping, Eq.(7)
if flags.do_mrs && flags.do_regular && kv.mrsmsp > 0
  
  angbelow = -90:kv.rangsamp:min(rang)-kv.rangsamp;
  angabove = max(rang)+kv.rangsamp:kv.rangsamp:265;
  rang = [angbelow,rang,angabove];
  ri = [zeros(length(angbelow),size(si,2)) ; si ; zeros(length(angabove),size(si,2))];

  mrs = kv.mrsmsp/cos(deg2rad(kv.lat)); % direction dependent scatter (derivation: const. length rel. to the circumferences of circles considered as cross sections of a unit sphere)

  x = 0:2*pi/length(rang):2*pi-2*pi/length(rang);
  kappa = 1/deg2rad(mrs)^2; % concentration parameter (~1/sigma^2 of normpdf)
  mrspdf = exp(kappa*cos(x)) / (2*pi*besseli(0,kappa)); % von Mises PDF 
  for tt = 1:size(ri,2)
    ri(:,tt) = pconv(ri(:,tt),mrspdf(:));
  end
    
else
  ri = si;
end
  

end


