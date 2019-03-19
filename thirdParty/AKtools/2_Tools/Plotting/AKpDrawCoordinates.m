function AKpDrawCoordinates(co, type, marker_size, line_width)
% AKpDrawCoordinates(co, type, marker_size, line_width)
% draws the coordinate axis in 3D plots.
% 
% See AKplotDemo.m for examples 
%
% I N P U T
% co          - strct with fields x, y,z (unit vectors pointing in positive
%               x, y, and z direction) and L (scalar that gives the axis
%               length)
% type        - 1: lines are drawn, marker indicates positive x,y and z
%               2: additional circles and lines mark intermediate angles
% marker_size - Size of markers that show positive x, y, and z direction
% line_width  - Line width of x, y, and z axis
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin,
% DFG research unit 'SEACEN', 7/2012

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

if ~exist('marker_size', 'var')
    % this is the size for the marker that indicates positive x,y and z
    marker_size = 10;
end

colors = AKcolors;

switch type
    case 1
        hold on
        % draw x axis and point for markin positive x
        plot3([-co.x(1)*co.L co.x(1)*co.L], ...
              [-co.x(2)*co.L co.x(2)*co.L], ...
              [-co.x(3)*co.L co.x(3)*co.L], 'color', colors.rgb(4,:), 'LineWidth', line_width)
        plot3(co.x(1)*co.L, co.x(2)*co.L, co.x(3)*co.L, '.', 'color', colors.rgb(4,:), 'Markersize', marker_size)
        % draw y axis and point for markin positive y
        plot3([-co.y(1)*co.L co.y(1)*co.L], ...
              [-co.y(2)*co.L co.y(2)*co.L], ...
              [-co.y(3)*co.L co.y(3)*co.L], 'color', colors.rgb(7,:), 'LineWidth', line_width)
        plot3(co.y(1)*co.L, co.y(2)*co.L, co.y(3)*co.L, '.', 'color', colors.rgb(7,:), 'Markersize', marker_size)
        % draw z axis and point for markin positive z
        plot3([-co.z(1)*co.L co.z(1)*co.L], ...
              [-co.z(2)*co.L co.z(2)*co.L], ...
              [-co.z(3)*co.L co.z(3)*co.L], 'color', colors.rgb(3,:), 'LineWidth', line_width)
        plot3(co.z(1)*co.L, co.z(2)*co.L, co.z(3)*co.L, '.', 'color', colors.rgb(3,:), 'Markersize', marker_size)
        hold off
    case 2
        hold on
        % draw circles
        x = linspace(0, 2*pi, 256);
        y = zeros(size(x));
        plot3(cos(x)*co.L, sin(x)*co.L, y, 'color', [.7 .7 .7], 'LineWidth', line_width)
        plot3(cos(x)*co.L, y, sin(x)*co.L, 'color', [.7 .7 .7], 'LineWidth', line_width)
        plot3(y, cos(x)*co.L, sin(x)*co.L, 'color', [.7 .7 .7], 'LineWidth', line_width)
        % draw angles
        for aa = [22.5 45 67.5 112.5 135 157.5]
            plot3(cosd(aa)*[co.L -co.L], sind(aa)*[co.L -co.L], [0 0], 'color', [.7 .7 .7], 'LineWidth', line_width)
            plot3(cosd(aa)*[co.L -co.L], [0 0], sind(aa)*[co.L -co.L], 'color', [.7 .7 .7], 'LineWidth', line_width)
            plot3([0 0], cosd(aa)*[co.L -co.L], sind(aa)*[co.L -co.L], 'color', [.7 .7 .7], 'LineWidth', line_width)
        end
        % draw x axis and point for markin positive x
        plot3([-co.x(1)*co.L co.x(1)*co.L], ...
              [-co.x(2)*co.L co.x(2)*co.L], ...
              [-co.x(3)*co.L co.x(3)*co.L], 'color', colors.rgb(4,:), 'LineWidth', line_width)
        plot3(co.x(1)*co.L, co.x(2)*co.L, co.x(3)*co.L, '.', 'color', colors.rgb(4,:), 'Markersize', marker_size)
        % draw y axis and point for markin positive y
        plot3([-co.y(1)*co.L co.y(1)*co.L], ...
              [-co.y(2)*co.L co.y(2)*co.L], ...
              [-co.y(3)*co.L co.y(3)*co.L], 'color', colors.rgb(7,:), 'LineWidth', line_width)
        plot3(co.y(1)*co.L, co.y(2)*co.L, co.y(3)*co.L, '.', 'color', colors.rgb(7,:), 'Markersize', marker_size)
        % draw z axis and point for markin positive z
        plot3([-co.z(1)*co.L co.z(1)*co.L], ...
              [-co.z(2)*co.L co.z(2)*co.L], ...
              [-co.z(3)*co.L co.z(3)*co.L], 'color', colors.rgb(3,:), 'LineWidth', line_width)
        plot3(co.z(1)*co.L, co.z(2)*co.L, co.z(3)*co.L, '.', 'color', colors.rgb(3,:), 'Markersize', marker_size)
        hold off
end
