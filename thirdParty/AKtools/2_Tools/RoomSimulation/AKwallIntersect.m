% intersect = AKwallIntersect(planeNormal,planePoint,lineA,lineB, L)
% cecks if there is an intersection between a line segment and a finite,
% rectangular plane
%
% for example use see AKism.m
%
% I N P U T
% planeNormal - normal vector of the plane [x y z]
% planePoint  - a point on the plane [x y z]
% lineA       - first end point of the line segment [x y z]
% lineB       - second end point of the line segment [x y z]
% L           - room dimensions according to AKism.m [x y z]
%
% O U T P U T
% intersect   - is the point of interection [x y z] if it exists, and false
%               otherwise
%
% 2017/05 - fabian.brinkmann@tu-berlin.de

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
function intersect = AKwallIntersect(planeNormal,planePoint,lineA,lineB, L)

% vector pointing in the line direction
lineDir = lineB-lineA;

% check if the line is not parallel to the wall
pN_dot_lD = dot(planeNormal,lineDir);

if abs(pN_dot_lD) > 1e-7
    
    % check if the line segment intersects the infinite plane
    pN_dot_pl = dot(planeNormal, planePoint-lineA);
    
    % get the intersection point
    d          = pN_dot_pl / pN_dot_lD;
    
    if d >= 0 && d <= 1
        intersect  = lineA+ d.*lineDir;
        
        % check if the intersection point is part of the wall
        if intersect(1)>=-1e-6 && intersect(1)<=L(1)+1e-6 && ...
                intersect(2)>=-1e-6 && intersect(2)<=L(2)+1e-6 && ...
                intersect(3)>=-1e-6 && intersect(3)<=L(3)+1e-6
            
            % clip intersect to the room dimensions
            intersect(1) = max(0, intersect(1));
            intersect(1) = min(L(1), intersect(1));
            intersect(2) = max(0, intersect(2));
            intersect(2) = min(L(2), intersect(2));
            intersect(3) = max(0, intersect(3));
            intersect(3) = min(L(3), intersect(3));
            
        else
            intersect = false;
        end
    else
        intersect = false;
    end
    
else
    intersect = false;
end