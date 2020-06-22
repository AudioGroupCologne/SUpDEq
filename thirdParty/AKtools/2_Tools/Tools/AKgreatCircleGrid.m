function [hrtf_grid_deg, act_ang_GCD, act_ang, weights] = AKgreatCircleGrid(el, max_ang, fit, do_plot, res_ang)
% [hrtf_grid_deg, act_ang_GCD, act_ang] = AKgreatCircleGrid(el, max_ang, fit, do_plot)
%
% Constructs a spherical grid with the following criterion:
% great arc distance between neighboring points of the same elevation has
% to be smaller or equal to max_ang degrees.
%
% AKgreatCircleGrid(-90:10:90, 10)
% creates a grid with 10 degree spacing in azimuth and elevation
%
%
% INPUT:
% el      - vector containing all elevation in degree
%           (row- or column-vector, default is 90:-2:-90) 
% max_ang - maximum great circle distance between to neighboring points of
%           the same elevation (default is 2)
% fit     - great circle distance is choosen to include an azimuth  point
%           each n*fit degrees. If fit = 90, azimuth values at 0, 90, 180,
%           and 270 degrees are included (i.e. median and frontal plane are
%           included). If fit = 0, or fit = 360 no specific spherical cross
%           section is included in the grid and the actual distances
%           between neighboring azimuth values are closest to the desired
%           great circle distance (default = 90).
% do_plot - plot resulting grid (default = 0)
% res_ang - angular resolution, e.g., if res_ang = 1 (default), the azimuth
%           values are multiples of 1.
%
% OUTPUT:
% hrtf_grid_deg - grid in degrees [azimut x elevation]
%                  azimuth   : 0:360 deg (0 in front of listener,
%                                         numeration counter clockwise)
%                  elevation : 90:-90 (0 in front of listener, 90 above)
% act_ang_GCD   - actual great circle distances used per elevation. Ordered
%                 according to 'el'
% act_ang       - actual angle used per elevation. Ordered according to
%                 'el'
% weights       - areas around the points in hrtf_grid_deg calculated from
%                 latitude-longitude-rectangles around them. weights are
%                 normalized to sum(weights) = 1.
%
%
% F. Brinkmann, Audio Communication Group, TU Berlin, 07/2013

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
if ~exist('el', 'var')
    el = 90:-2:-90;
end
if ~exist('max_ang', 'var')
    max_ang = 2;
end
if ~exist('fit', 'var')
    fit = 90;
end
if ~exist('do_plot', 'var')
    do_plot = 0;
end
if ~exist('res_ang', 'var')
    res_ang = 1;
end

% check input format
el = reshape(el, [numel(el) 1]);

if fit == 0
    fit = 360;
end
if rem(1, res_ang)
     warning('AKgreatCircleGrid:Input', 'results can contain errors if 1/res_ang is NOT an integer numer')
end

% calculate delta phi to meet the criterion
% (according to Bovbjerg et al. 2000: Measuring the head related transfer
% functions of an artificial head with a high directional resolution,
% R.W. Sinnott, "Virtues of the Haversine", Sky and Telescope, vol. 68, no.
% 2, 1984, p. 159)
d_phi = 2*asind( sind(max_ang/2) ./ cosd(el) );
% correct values at the poles
d_phi(abs(el)==90) = 360;
% round to desired angular resolution
d_phi = floor( d_phi / res_ang ) * res_ang;

% adjust delta phi to assure equal spacing on sphere -> mod(360, d_phi)!=0
% or quarter sphere -> mod(90, d_phi)!=0
% (if equally spaced on on quarter sphere (fit = 90), median, horizontal
% and frontal plane are included in measurements).
% this operation is easier in degrees than in radians...
for n = 1:length(d_phi)
    if abs(el(n)) ~= 90
        while mod(fit, d_phi(n))
            d_phi(n) = round( (d_phi(n) - res_ang) / res_ang) * res_ang;
        end
    else
        % irregularity at north and south pole
        d_phi(n) = 360;
    end
end
act_ang = d_phi;
clearvars n

% calculate great circle angle that is actually used in the grid
% (R.W. Sinnott, "Virtues of the Haversine", Sky and Telescope, vol. 68, no. 2, 1984, p. 159)
act_ang_GCD = 2*asind(sqrt(cosd(el).^2.*sind(d_phi/2).^2));


% construct pre-grid
hrtf_grid = [];
m = 1;
for n = 1:length(d_phi)
    tmp = 0:d_phi(n):360-d_phi(n);
    hrtf_grid(m:m+length(tmp)-1, 1) = tmp;
    hrtf_grid(m:m+length(tmp)-1, 2) = el(n);
    m = m + length(tmp);
end

% final grid in degree
hrtf_grid_deg = hrtf_grid;

% estimated area weights using lat-long rectangles
weights       = nan(size(hrtf_grid_deg, 1), 1);
[el_sort, id] = sort(el);
act_ang_sort  = act_ang(id);

for n = 1:length(act_ang)
    
    if length(act_ang) == 1
        weight = 1;
    elseif el_sort(n) == -90
        el_range = [-90 mean(el_sort(n:n+1))];
        weight   = area_quad(el_range(1), -180, el_range(2), 180);
    elseif el_sort(n) == 90
        el_range = [90 mean(el_sort(n-1:n))];
        weight   = area_quad(el_range(1), -180, el_range(2), 180);
    else
        
        if n == 1
            el_diff  = ( el_sort(n+1)-el_sort(n) ) / 2;
            el_range = el_sort(n) + [-el_diff el_diff];
            weight   = area_quad(el_range(1), 0, el_range(2), act_ang_sort(n));
        elseif n == length(act_ang)
            el_diff  = ( el_sort(n)-el_sort(n-1) ) / 2;
            el_range = el_sort(n) + [-el_diff el_diff];
            weight   = area_quad(el_range(1), 0, el_range(2), act_ang_sort(n));
        else
            el_range = [mean(el_sort(n-1:n)) mean(el_sort(n:n+1))];
            weight   = area_quad(el_range(1), 0, el_range(2), act_ang_sort(n));
        end
        
    end
    
    weights( hrtf_grid_deg(:,2) == el_sort(n) ) = weight;
    
end

weights = weights / sum(weights);

% plot
if do_plot
    % convert to grid that suits matlabs plotting routine
    hrtf_grid(:, 2) = 90-hrtf_grid(:, 2);
    hrtf_grid       = degtorad(hrtf_grid);
    
    % plot final grid
    fWidth  = 20;
    fHeight = 10;
    hFigureHandle = figure;
    set(hFigureHandle,'PaperUnits', 'centimeters');
    set(hFigureHandle,'Units', 'centimeters');
    set(hFigureHandle, 'PaperSize', [fWidth fHeight]);
    set(hFigureHandle,'PaperPosition', [.1 .1 fWidth-.1 fHeight-.1]);
    set(hFigureHandle,'Position', [0 0 fWidth fHeight]);
    set(hFigureHandle, 'color', [1 1 1])
    
    subplot(1,2,1)
    tmp = find(hrtf_grid(:, 1) <= pi);
    [hrtf_x, hrtf_y, hrtf_z] = sph2cart(hrtf_grid(tmp, 1), hrtf_grid(tmp, 2)+pi/2, ones(length(tmp), 1));
    
    plot3(hrtf_x, hrtf_y, hrtf_z, '.k', 'MarkerSize', 2)
    axis square
    axis([-1.1 1.1 -1.1 1.1 -1.1 1.1])
    set(gca, 'xTick', [-1 0 1], 'yTick', [-1 0 1], 'zTick', [-1 0 1])
    view([0 0])
    title ('Side view')
    
    subplot(1,2,2)
    tmp = find(hrtf_grid(:, 2) <= pi/2);
    [hrtf_x, hrtf_y, hrtf_z] = sph2cart(hrtf_grid(tmp, 1), hrtf_grid(tmp, 2)+pi/2, ones(length(tmp), 1));
    
    plot3(hrtf_x, hrtf_y, hrtf_z, '.k', 'MarkerSize', 2)
    axis square
    axis([-1.1 1.1 -1.1 1.1 -1.1 1.1])
    set(gca, 'xTick', [-1 0 1], 'yTick', [-1 0 1], 'zTick', [-1 0 1])
    view([0 90])
    title ('Top view')
end

clearvars tmp hrtf_x hrtf_y hrtf_z hFigureHandle fWidth fHeight m

end

% function to replace Matlab's areaquad from the Mapping Toolbox
function area = area_quad(lat1, lon1, lat2, lon2)
    area = 2*pi * abs( sind(lat1)-sind(lat2) ) * abs( lon1-lon2 ) / 360;
end
