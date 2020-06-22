% Demonstrates the usage of AKsubGrid which can be called to extract
% certain points from a spherical sampling grid
%
% 03/2016 - fabian.brinkmann@tu-berlin.de

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
close all; clear; clc

% generate a spherical sampling grid
gridIn = AKgreatCircleGrid(90:-5:-90, 5, 90, 0);

% suppress unused variables warning for demo script
%#ok<*ASGLU>

%% select any point from the grid, we use a point in frontal direction and
%  additionally select all points within a great circle distance of 15
%  degree
[id, subGrid, subGridHor] = AKsubGrid(gridIn, 'any', [0 0], 15, true);


%% select points on a transverse plane and plot
[id, subGrid, subGridHor] = AKsubGrid(gridIn, 'transverse', 30, 0, true);


%% select points on a sagittal plane and plot
[id, subGrid, subGridHor] = AKsubGrid(gridIn, 'sagittal', 30, 5, true);

%% select points on a corconal plane and plot
[id, subGrid, subGridHor] = AKsubGrid(gridIn, 'corconal', 45, 5, true);

%% generate a horizontal plane
[id, subGrid, subGridHor] = AKsubGrid(2, 'transverse', 0, 0, true);

%% generate a median plane
[id, subGrid, subGridHor] = AKsubGrid(2, 'sagittal', 0, 0, true);

%% generate a frontal plane
[id, subGrid, subGridHor] = AKsubGrid(2, 'corconal', 0, 0, true);

%% generate a Gauss-like grid with in vertical polar coordinates
[id, subGrid, subGridHor] = AKsubGrid(10, 'transverse', -90:10:90, 0, true);

%% generate a Gauss-like grid with in horizontal polar coordinates
[id, subGrid, subGridHor] = AKsubGrid(10, 'sagittal', -90:10:90, 0, true);

%% generate a Gauss-like grid with in frontal polar coordinates
[id, subGrid, subGridHor] = AKsubGrid(10, 'corconal', -90:10:90, 0, true);
