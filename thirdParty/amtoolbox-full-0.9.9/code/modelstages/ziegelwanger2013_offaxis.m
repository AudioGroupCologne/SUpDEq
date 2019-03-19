function y=ziegelwanger2013_offaxis(p,x)
%ziegelwanger2013_offaxis XXX
%   Usage: y=ziegelwanger2013_offaxis(p,x) 
%
%   Input parameters:
%       p: off-axis time-of-arrival model parameters [SI-units]
%       x: HRTF direction (azimuth,elevation) [rad]
%   Output parameters:
%       y: time-of-arrival [s]
%
%   toa=ZIEGELWANGER2013_OFFAXIS(p,x) calculates time-of-arrivals for given
%   model parameters (p) and directions (x) with an off-axis time-of-arrival
%   model.
%
%   See also: ziegelwanger2013, ziegelwanger2013_onaxis,
%   data_ziegelwanger2013, exp_ziegelwanger2013
%
%   References:
%     P. Majdak and H. Ziegelwanger. Continuous-direction model of the
%     broadband time-of-arrival in the head-related transfer functions. In
%     ICA 2013 Montreal, volume 19, page 050016, Montreal, Canada, 2013. ASA.
%     
%     H. Ziegelwanger and P. Majdak. Modeling the broadband time-of-arrival
%     of the head-related transfer functions for binaural audio. In
%     Proceedings of the 134th Convention of the Audio Engineering Society,
%     page 7, Rome, 2013.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/ziegelwanger2013_offaxis.php

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

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria

r=p(1); %............. sphere radius [m]
xM=p(2); %............ x-coordinate of the sphere center [m]
yM=p(3); %............ y-coordinate of the sphere center [m]
zM=p(4); %............ z-coordinate of the sphere center [m]
delay=p(5); %......... constant dely [s]
phi_ear=p(6); %....... position of the ear (azimuth angle) [rad]
theta_ear=p(7); %..... position of the ear (elevation angle) [rad]

M=sqrt(xM^2+yM^2+zM^2);

beta=acos(-cos(x(:,2)).*(xM*cos(x(:,1))+yM*sin(x(:,1)))-zM*sin(x(:,2)));
s2=-r+M*cos(beta)+sqrt(r^2+M^2*cos(beta).^2+2*M*r);
gamma=pi-beta-acos((2*M^2+2*M*r-2*r*s2-s2.^2)/(2*M^2+2*M*r));
if M==0
    s1=zeros(size(x,1),1);
else
    s1=M*cos(beta)./(2*(M+r).*tan(gamma/2));
end

y=1/343*((r* ...
    ((sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
    (1-sin(theta_ear).*sin(x(:,2))-cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))+ ...
    (-sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
    (1+acos(sin(theta_ear).*sin(x(:,2))+cos(theta_ear)*cos(x(:,2)).*cos(phi_ear-x(:,1)))-pi/2))) ...
    +s1+s2) ...
    +delay-(M+r)/343;

end
