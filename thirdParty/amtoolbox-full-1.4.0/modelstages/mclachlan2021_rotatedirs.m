function [dirs_rot] = mclachlan2021_rotatedirs(dirs,theta,type)
%MCLACHLAN2021_ROTATEDIRS Rotate a set of coordinates on a sphere to a new set of coordinates
%   Usage: [dirs_rot] = rotateDirs(dirs,theta,type);
%
%   Input parameters:
%       dirs       : original directions in cartesian coordinates [x,y,z;...]
%       theta      : rotation amount (in degrees)
%       type       : type of rotation, options are 'yaw','pitch','roll'
%
%   Output parameters:
%       dirs_rot   : rotated directions in cartesian coordinates [x,y,z;...]
%
%   See also: mclachlan2021, mclachlan2021_preproc
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/mclachlan2021_rotatedirs.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal M-Image
%   #Author: Glen McLachlan (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% convert to rads, head rotation = -source rotation
theta = -theta/180*pi;

% rotate coordinates according to type of rotation
if strcmp(type,'yaw')
    dirs_rot(:,1) = cos(theta).*dirs(:,1) - sin(theta).*dirs(:,2) ;
    dirs_rot(:,2) = sin(theta).*dirs(:,1) + cos(theta).*dirs(:,2) ;
    dirs_rot(:,3) = dirs(:,3) ;
elseif strcmp(type,'pitch')
    dirs_rot(:,1) = cos(theta).*dirs(:,1) + sin(theta).*dirs(:,3) ;
    dirs_rot(:,2) = dirs(:,2);
    dirs_rot(:,3) = -sin(theta).*dirs(:,1) + cos(theta).*dirs(:,3) ;
elseif strcmp(type,'roll')
    dirs_rot(:,1) = dirs(:,1);
    dirs_rot(:,2) = cos(theta).*dirs(:,2) - sin(theta).*dirs(:,3);
    dirs_rot(:,3) = sin(theta).*dirs(:,2) + cos(theta).*dirs(:,3);

else
    error('Unrecognised rotation type');
end


