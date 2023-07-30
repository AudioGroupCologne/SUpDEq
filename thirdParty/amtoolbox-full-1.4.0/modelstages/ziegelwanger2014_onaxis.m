function y=ziegelwanger2014_onaxis(p,x)
%ZIEGELWANGER2014_ONAXIS   On-axis time-of-arrival model
%   Usage: y=ziegelwanger2014_onaxis(p,x)
%
%   Input parameters:
%       p: on-axis model parameters [SI-units]
%       x: HRTF direction (azimuth,elevation) [rad]
%   Output parameters:
%       y: time-of-arrival [s]
%
%   toa=ZIEGELWANGER2014_ONAXIS(p,x) calculates time-of-arrivals (TOAs) for
%   given model parameters (p) and directions (x) with an on-axis
%   time-of-arrival model.
%
%   See also: ziegelwanger2014, ziegelwanger2014_offaxis,
%   data_ziegelwanger2014, exp_ziegelwanger2014
%
%   References:
%     H. Ziegelwanger and P. Majdak. Modeling the direction-continuous
%     time-of-arrival in head-related transfer functions. J. Acoust. Soc.
%     Am., 135:1278--1293, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/ziegelwanger2014_onaxis.php


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
phi_ear=p(2); %....... position of the ear (azimuth angle) [rad]
theta_ear=p(3); %..... position of the ear (elevation angle) [rad]
delay=p(4); %......... constant delay [s]

y=r/340.*( ...
       (sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
       (1-sin(theta_ear).*sin(x(:,2))-cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))+ ...
       (-sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
       (1+acos(sin(theta_ear).*sin(x(:,2))+cos(theta_ear)*cos(x(:,2)).*cos(phi_ear-x(:,1)))-pi/2))+ ...
       delay-r/340;
end


