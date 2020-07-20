%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function d = supdeq_haversine(AzEl1, AzEl2, radius) 
%
% This function calculates the great circle distance between two points on
% a sphere by applying the inverse haversine function 
%
% Output:
% d                     - Great circle distance d between two points on the sphere
%
% Input:
% AzEl1 / AzEl2         - Spatial sampling point (e.g., [0,90]), where the first 
%                       value is the azimuth and the second the elevation (both in degree).
%                       Azimuth in degree (0=front, 90=left, 180=back, 270=right)
%                       (0 points to positive x-axis, 90 to positive y-axis)
%                       Elevations in degree (0=North Pole, 90=front, 180=South Pole)
%                       (0 points to positive z-axis, 180 to negative z-axis)
% radius                - Head/Sphere radius in m
%                       Default: 0.0875 m
%
% Dependencies: -
%
% (C) 2020 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function d = supdeq_haversine(AzEl1, AzEl2, radius) 

if nargin < 3 || isempty(radius)
    radius = 0.0875;
end

%% Transform sampling points to radiant in matlab coordinates system

AzEl1 = AzEl1*pi/180;
AzEl1(:,2) = pi/2 - AzEl1(:,2);

AzEl2 = AzEl2*pi/180;
AzEl2(:,2) = pi/2 - AzEl2(:,2);

%% Apply formula

%Get distance between angular values
dAz = AzEl2(:,1) - AzEl1(:,1);
dEl = AzEl2(:,2) - AzEl1(:,2);

%Get great circular distance
%d = 2 * radius * asin( sqrt( sin(dAz/2).^2 + cos(AzEl1(:,1)) .* cos(AzEl2(:,1)) .* sin(dEl/2).^2));
d = 2 * radius * asin( sqrt( sin(dEl/2).^2 + cos(AzEl1(:,2)) .* cos(AzEl2(:,2)) .* sin(dAz/2).^2)); %According to 3D Tune In

end

