function y=ziegelwanger2014_offaxis(p,x)
%ZIEGELWANGER2014_OFFAXIS   Off-axis time-of-arrival model
%   Usage: y=ziegelwanger2014_offaxis(p,x) 
%
%   Input parameters:
%       p: off-axis time-of-arrival model parameters [SI-units]
%       x: HRTF direction (azimuth,elevation) [rad]
%   Output parameters:
%       y: time-of-arrival [s]
%
%   toa=ZIEGELWANGER2014_OFFAXIS(p,x) calculates time-of-arrivals for given
%   model parameters (p) and directions (x) with an off-axis time-of-arrival
%   model.
%
%   See also: ziegelwanger2014, ziegelwanger2014_onaxis,
%   data_ziegelwanger2014, exp_ziegelwanger2014
%
%   References:
%     H. Ziegelwanger and P. Majdak. Modeling the direction-continuous
%     time-of-arrival in head-related transfer functions. J. Acoust. Soc.
%     Am., 135:1278--1293, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/ziegelwanger2014_offaxis.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA M-Signal M-Optimization
%   #Author: Harald Ziegelwanger (2014), Acoustics Research Institute, Vienna, Austria
%   #Author: Robert Baumgartner (2018)
%   #Author: Laurin Steidle (2018)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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

y=1/340*((r* ...
    ((sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
    (1-sin(theta_ear).*sin(x(:,2))-cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))+ ...
    (-sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
    (1+acos(sin(theta_ear).*sin(x(:,2))+cos(theta_ear)*cos(x(:,2)).*cos(phi_ear-x(:,1)))-pi/2))) ...
    +s1+s2) ...
    +delay-(M+r)/340;

end


