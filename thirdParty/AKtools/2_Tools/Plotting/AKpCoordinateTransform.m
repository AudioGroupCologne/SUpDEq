function co = AKpCoordinateTransform(az, el, type)
% function co = hp_coordinate_transformation(az, el, type)
%
% transform different coordinate conventions to matlab default
% if az containes values bigger than 2pi, it is assumed that the unit of az
% is degree. In this case it is converted to radians
%
% See AKplotDemo.m for examples
%
% type
% 1: Matlab default
% 0 <= az < 360 (0=source in front 90 = source to the left of listener)
% 90 <= el <= -90 (90=north pole, 0=front, -90=south pole)
% positive x is at (0|90) (az|el)
% positive y is at (90|90)
% positive z is at (0|0);
%
% 2: Mathematical
% 0 <= az < 360 (0=source in front 90 = source to the left of listener)
% 0 <= el <= 180 (0=north pole, 90=front, 180=south pole)
% positive x is at (0|90) (az|el)
% positive y is at (90|90)
% positive z is at (0|0);

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

if max(max(az > 2*pi))
    az = deg2rad(az);
    if ~islogical(el)
        el = deg2rad(el);
    else
        el = 0;
    end
end

switch type
    case 1 % matlab default, see doc sph2cart
        % coordinate transformation
        co.az = az;
        co.el = el;
        % definition of axis normals for balloon plots
        co.x = [1 0 0]';
        co.y = [0 1 0]';
        co.z = [0 0 1]';
        % x and y axis for spherical planar plots
        co.xtick      = 1:45:361;
        co.xticklabel = 0:45:360;
        co.ytick      = 1:45:181;
        co.yticklabel = [90:-45:0 -45 -90];
    case 2 % mathematical correct
        
        % coordinate transformation
        co.az = az;
        id1     = el > pi/2;
        id2     = el <= pi/2;
        el(id1) = -(el(id1) - pi/2);
        el(id2) = abs(el(id2) - pi/2);
        co.el   = el;
        % definition of axis normals for balloon plots
        co.x = [1 0 0]';
        co.y = [0 1 0]';
        co.z = [0 0 1]';
        % x and y axis for spherical planar plots
        co.xtick      = 1:45:361;
        co.xticklabel = 0:45:360;
        co.ytick      = 1:45:181;
        co.yticklabel = 0:45:180;
end
