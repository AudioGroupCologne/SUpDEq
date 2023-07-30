function [lat,pol]=sph2horpolar(azi,ele)
%   SPH2HORPOLAR from spherical to horizontal-polar coordinate system
%
%   Usage:    [lat,pol] = sph2horpolar(azi,ele)
%
%   Input parameters:
%     azi     : azimuth in deg
%     ele     : elevation in deg
%
%   Output parameters:
%     lat     : lateral angle in deg, [-90 deg, +90 deg]
%     pol     : polar angle in deg, [-90 deg, 270 deg]
%
%   SPH2HORPOLAR(...) converts spherical coordinates (azimuth and 
%   elevation angles) into the coordinates of the horizontal-polar system 
%   (lateral and polar angles).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/sph2horpolar.php


%   #Author: Peter L. Sondergaard

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

azi=mod(azi+360,360);
ele=mod(ele+360,360);

razi = deg2rad(azi);
rele = deg2rad(ele);
rlat=asin(sin(razi).*cos(rele));
rpol=real(asin(sin(rele)./cos(rlat)));
rpol(cos(rlat)==0)=0;
pol = rad2deg(rpol);
lat = rad2deg(rlat);

idx = find(razi>pi/2 & razi < 3*pi/2 & (rele < pi/2 | rele > 3*pi/2));
pol(idx)=180-pol(idx);
idx = find(~(razi>pi/2 & razi < 3*pi/2) & rele > pi/2 & rele < 3*pi/2);
pol(idx)=180-pol(idx);

end


