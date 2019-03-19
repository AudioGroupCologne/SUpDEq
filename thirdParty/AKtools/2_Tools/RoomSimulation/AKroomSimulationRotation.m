% [azRot, elRot] = AKroomSimulationRotation(az, el, rotAz, rotEl)
% applies a rotation given by azRot and elRot to incident and exit angels
% at receiver and sources given by az and el.
% The source/reciever rotation is defined about the global z-axis
% (left/right rotation), and the local y-axis of the source/receiver
% (up/down rotation). This way the position after rotation does not depend
% on the order in which the rotations are done. The rotated values azRot
% and elRot are obtained by the inverse rotation being applied to the input
% values az and el:
% The source/receiver rotation is given by Rz*Ry*X, where Raz and Rel are
% rotation matrixes and X is a point on the source/receiver in Cartesian
% coordinates. The inverse rotation is given by (Rz*Ry)'*Y, where Y is the
% direction of the incoming/outgoiung sound path expressed in Cartesian
% coordinates.
%
% to be called from AKroomSimulation - see this function for use cases
%
% I N P U T
% az    - azimuth values of incoming/outgoing sound paths [N x 1]
% el    - elevation values of incoming/outgoing sound paths [N x 1]
% rotAz - azimuth rotation values [1 x M]
% rotEl - elevation rotation values [1 x M]
%
% O U T P U T
% azRot - N by M matrix of rotated azimut values
% elRot - N by M matrix of rotated elevation values
%
% 2018/01 - fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [azRot, elRot] = AKroomSimulationRotation6(az, el, rotAz, rotEl)

% check input format
try
    az    = reshape(az,    [numel(az)    1]);
    el    = reshape(el,    [numel(az)    1]);
    rotAz = reshape(rotAz, [1 numel(rotAz)]);
    rotEl = reshape(rotEl, [1 numel(rotAz)]);
catch
    error('AKroomSimulationRotation:Input', 'Input data has wrong size: ''az'' and ''el'' must be [N x 1], ''rotAz'' and ''rotEl'' must be [1 x M].')
end

% convert input points to Cartesian coordinates
[x, y, z] = sph2cart(az/180*pi, el/180*pi, ones(size(az)));
X         = [x'; y'; z'];

% apply the inverse rotation
azRot = nan(size(az,1), size(rotAz,2));
elRot = azRot;

for nn = 1:numel(rotAz)
    % rotation matrix about z-Axis
    rAz = [cosd(rotAz(nn)) -sind(rotAz(nn)) 0;...
           sind(rotAz(nn))  cosd(rotAz(nn))  0;...
           0                0                1];

    % rotation matrix about y-Axis
    rEl = [cosd(-rotEl(nn)) 0 sind(-rotEl(nn));...
           0                 1 0;...
           -sind(-rotEl(nn)) 0 cosd(-rotEl(nn))];
    
    % apply inverse rotation to the incoming/outgoing sound paths
    % (rAz^-1 = rAz', and (rAz*rEl)' = rEl'*rAz')
    Xr = (rAz*rEl)' * X;

    % convert to spherical coordinates
    [azRot(:,nn), elRot(:,nn)] = cart2sph(Xr(1,:)', Xr(2,:)', Xr(3,:)');
end

% rad to deg
azRot = azRot/pi*180;
azRot = mod(azRot, 360);
elRot = elRot/pi*180;