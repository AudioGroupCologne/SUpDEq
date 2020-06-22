function [X, Y, Z, data] = AKpDataInterpolation(co, data)
% interpolates spherical data with an arbitrary grid to a full sphere with
% 1 degree resolution. Works only with the matlab default definition of
% coordinates (see doc sph2cart). Warning: Simple interpolation is
% performed that might distort the data!
%
% Only to be called from AKp.m
% See AKplotDemo.m for examples
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
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License.

% transpose az and el to obtain row vectors
if size(co.az,1) > size(co.az,2)
    co.az = co.az';
end
if size(co.el,1) > size(co.el,2)
    co.el = co.el';
end

% show original and interpolated data (for debugging only)
do_plot = 0;

if do_plot
    scatter3(co.az, co.el, data) %#ok<*UNRCH>
end

% the interpolation is done one 3D plane, where az is x axis, el is y axis
% and data is z - axis. For propper interpolation close to the edges of the
% plane, duplicate version of the original data has to be moved to the left
% and right and mirrored above and below. This is done with the full data,
% to make sure that it works in any case. It therefore might be a bit
% slower. This means we need to add 8 surrounding versions of the original
% data.

% left and right
az_i   = [co.az co.az-2*pi co.az+2*pi];
el_i   = [co.el co.el co.el];
data_i = [data data data];

% mirror everthing to the top
id1 = el_i >= 0;
id2 = el_i < 0;

el1 = el_i;
el1(id1) = el1(id1) + 2*(pi/2-el1(id1));
el1(id2) = abs(el1(id2)) + pi;

% mirror everything to the botton
el2 = el_i;
el2(id1) = -el2(id1) - pi;
el2(id2) = el2(id2) -2*(pi/2+el2(id2));

% insert mirrored versions
az_i   = [az_i az_i az_i];
el_i   = [el_i el1 el2];
data_i = [data_i data_i data_i];

% get interpolator object
warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
data_interpolator  = scatteredInterpolant(az_i',el_i',data_i');
warning('on', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')

% get grid for interpolation
az = (0:360)*pi/180;
el = (90:-1:-90)*pi/180;
[az, el] = meshgrid(az, el);

% get interpolated data
data = data_interpolator(az, el);

% transform to cartesian coordinates
[X, Y, Z] = sph2cart(az, el, data);

if do_plot
    hold on
    mesh(az, el, data)
    hold off
end
