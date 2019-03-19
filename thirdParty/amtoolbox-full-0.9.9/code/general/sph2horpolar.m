function [lat,pol]=sph2horpolar(azi,ele)
% SPH2HORPOLAR from spherical to horizontal-polar coordinate system
%
% Usage:    [lat,pol] = sph2horpolar(azi,ele)
%
% Input parameters:
%     azi     : azimuth in deg
%     ele     : elevation in deg
%
% Output parameters:
%     lat     : lateral angle in deg, [-90°,+90°]
%     pol     : polar angle in deg, [-90°,270°]
%
%   geo2horpolar(...) converts spherical coordinates (azimuth and 
%   elevation) into the coordinates of the horizontal-polar system, i.e.,
%   lateral and polar angles.
%
% AUTHOR: Piotr Majdak, Robert Baumgartner
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/sph2horpolar.php

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

azi=mod(azi+360,360);
ele=mod(ele+360,360);

razi = deg2rad(azi);
rele = deg2rad(ele);
rlat=asin(sin(razi).*cos(rele));
rpol=asin(sin(rele)./cos(rlat));
rpol(cos(rlat)==0)=0;
pol = rad2deg(rpol);
lat = rad2deg(rlat);

idx = find(razi>pi/2 & razi < 3*pi/2 & (rele < pi/2 | rele > 3*pi/2));
pol(idx)=180-pol(idx);
idx = find(~(razi>pi/2 & razi < 3*pi/2) & rele > pi/2 & rele < 3*pi/2);
pol(idx)=180-pol(idx);

end
