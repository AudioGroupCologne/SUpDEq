% [id, subGrid, subGridHor] = AKsubGrid(gridIn, subGrid, angle, tol, doPlot)
% extracts arbitrary points from a spatial spherical sampling grid
% extracts transverse, sagittal, and corcoanal planes
% generates transverse, sagittal, and corcoanal planes
%
% AKsubGrid(gridIn, 'any', [90 0])      -> extracts point to the right
% AKsubGrid(gridIn, 'transverse', 0)    -> extracts horizontal plane
% AKsubGrid(gridIn, 'sagittal', 0)      -> extracts median plane
% AKsubGrid(gridIn, 'corconal', 90)     -> extracts frontal plane
% AKsubGrid(10, 'transverse, 0)         -> generates horizontal plane
%                                          (10 degree resolution)
% AKsubGrid(10, 'transverse, -90:10:90) -> generates gauss like grid
%                                          (10 degree resolution)
%
% See AKsubGridDemo.m for more examples
%
%
% INPUT
% gridIn  - To extract points from a sampling grid, pass a spatial sampling
%           grid given by azimuth (first column) and elevation (second).
%           To generate a plane, pass a scalar that specifies the
%           resolution in degree.
%           Azimuth angles in degree ranging from 0 to 360 (0=front,
%           90=back). elevation angles in degree ranging from -90 to 90
%           (-90=below, 90=above)
% subGrid - 'any' will extract the point that is closed to angle (see
%                 below)
%           'transverse' will extract/generate points with same elevation,
%                        e.g. the orizontal plane.
%           'sagittal' will extract/generate points with same azimuth, e.g.
%                      the median plane, after transforming az and el into
%                      interaural polar coordinates. Here the azimuth
%                      ranges from -90 to 90 degrees and the elevation from
%                      -90 to 270 degrees.
%           'corconal' will extract/generate points with same absolute
%                      azimuth, e.g. the frontal plane with azimuths of 90
%                      and -90 (270) degree.
% angle   - specifies the desired angle [azimuth elevation] to extract if
%           subGrid='any'. In this case angles are specified in the format
%           given by gridIn. Otherwise this specifies the elevation angle
%           of a transverse plane, the azimuth angle of a sagittal plane, 
%           or the azimuth angle of a corconal plane respectively [degree]
% tol     - All grid points within the search range [angle-tol angle+tol]
%           are extracted (default=0.01 degree, only needed if points are
%           extracted from a grid)
% doPlot  - Plots extracted, and not extracted points (default=false)
%
%
% OUTPUT
% id         - indices of selectec points 
% subGrid    - points in the subGrid specifed by azimuth and elevation
% subGridHor - points in the subGrid given in interaural polar coordinates
%
%
% fabian.brinkmann@mailbox.org, Audio Communicatin Group, TU Berlin
% 01/2016

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
function [id, subGrid, subGridHor] = AKsubGrid(gridIn, subGrid, angle, tol, doPlot)

if ~exist('angle', 'var')
    if strcmpi(subGrid, 'any')
        angle = [90 0];
    else
        angle = 0;
    end
end
if ~exist('tol', 'var')
    tol = .01;
end
if ~exist('doPlot', 'var')
    doPlot = false;
end

if strcmpi(subGrid, 'any')
    
    if size(angle,1)~=1 || size(angle,2)~=2
        error('angle must be a 1x2 vector, i.e. [90 0]')
    end
    
    % find closest grid point to desired angle
    % (latidude = elevation; longitude = azimuth)
    if exist('distance.m', 'file')
        dMin = distance(angle(2), angle(1), gridIn(:,2), gridIn(:,1));
    else
        dMin = 2*asind( sqrt( sind( abs(angle(2)-gridIn(:,2)) ./ 2 ).^2 + cosd(angle(2)) * cosd(gridIn(:,2)) .* sind( abs(angle(1)-gridIn(:,1)) ./ 2 ).^2 ) );
    end
                
    [~,idMin] = min(dMin);
    
    angle = gridIn(idMin,:);
    
    % find all points within tolerance
    if exist('distance.m', 'file')
        d = distance(angle(2), angle(1), gridIn(:,2), gridIn(:,1));
    else
        d = 2*asind( sqrt( sind( abs(angle(2)-gridIn(:,2)) ./ 2 ).^2 + cosd(angle(2)) * cosd(gridIn(:,2)) .* sind( abs(angle(1)-gridIn(:,1)) ./ 2 ).^2 ) );
    end
            
    id      = find(d<=tol);
    subGrid = [gridIn(id,1) gridIn(id,2)];
    
    % transform to inteaural polar
    if exist('sph2hor', 'file')
        [azHor, elHor] = sph2hor(subGrid(:,1), subGrid(:,2));
        subGridHor     = [azHor elHor];
    else
        error('This needs the function ''sph2hor'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
    end
    
    
elseif strcmpi(subGrid, 'transverse')
    
    % extract points from grid
    if numel(gridIn) > 1
        % get minimum and maximum elevations
        elMax = min( 90, angle + abs(tol));
        elMin = max(-90, angle - abs(tol));
        
        % find desired grid points
        id = find(gridIn(:,2)<=elMax & gridIn(:,2)>=elMin);
        subGrid      = gridIn(id,1);
        subGrid(:,2) = gridIn(id,2);
        
        % sort
        [~, idSort] = sort(subGrid(:,1));
        id = id(idSort);
        subGrid(:,1) = subGrid(idSort,1);
        subGrid(:,2) = subGrid(idSort,2);
        
        % transform to inteaural polar
        if exist('sph2hor', 'file')
            [azHor, elHor] = sph2hor(subGrid(:,1), subGrid(:,2));
            subGridHor     = [azHor elHor];
        else
            error('This needs the function ''sph2hor'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
        end
        
    % generate grid
    else
        id      = [];
        
        % get azimuth
        az      = (0:gridIn:360)';
        if az(end) == 360
            az = az(1:end-1);
        end
        
        % construc the grid
        angle   = reshape( angle, numel(angle), 1 );
        subGrid = [repmat( az, numel(angle), 1 ) repelem( angle, numel(az), 1 )];
        
        % transform to inteaural polar
        if exist('sph2hor', 'file')
            [azHor, elHor] = sph2hor(subGrid(:,1), subGrid(:,2));
            subGridHor     = [azHor elHor];
        else
            error('This needs the function ''sph2hor'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
        end
    end
    
    
elseif strcmpi(subGrid, 'sagittal')
    
    % extract points from grid
    if numel(gridIn) > 1
        
        % transform coordinates to interaural polar
        if exist('sph2hor', 'file')
            [azHor, elHor] = sph2hor(gridIn(:,1), gridIn(:,2));
        else
            error('This needs the function ''sph2hor'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
        end
        
        % get minimum and maximum azimuths
        azMax = min( 90, angle+abs(tol));
        azMin = max(-90, angle-abs(tol));
        
        % find desired grid points
        id = find(azHor<=azMax & azHor>=azMin);
        subGridHor = [azHor(id)    elHor(id)];
        subGrid    = [gridIn(id,1) gridIn(id,2)];
        
        % sort
        [~, idSort] = sort(subGridHor(:,2));
        id    = id(idSort);
        subGridHor = [subGridHor(idSort,1) subGridHor(idSort,2)];
        subGrid    = [subGrid(idSort,1)    subGrid(idSort,2)];
   
    % generate grid
    else
        id         = [];
        
        % get elevation
        el         = (-90:gridIn:270)';
        if el(end) == 270
            el = el(1:end-1);
        end
        
        % construc the grid
        angle   = reshape( angle, numel(angle), 1 );
        subGridHor = [repelem( angle, numel(el), 1 ) repmat( el, numel(angle), 1 )];
        
        % transform to vertical polar
        if exist('hor2sph', 'file')
            [az, el] = hor2sph(subGridHor(:,1), subGridHor(:,2));
            subGrid     = [az el];
        else
            error('This needs the function ''hor2sph'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
        end 
    end
    
elseif strcmpi(subGrid, 'corconal')
    
    % extract points from grid
    if numel(gridIn) > 1
        % rotate, and transform coordinates to interaural polar
        if exist('sph2hor', 'file')
            [azHor, elHor] = sph2hor(mod(gridIn(:,1)+90,360), gridIn(:,2));
        else
            error('This needs the function ''sph2hor'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
        end
        
        % convert azimuth to original range
        azHor = 90-azHor;
        
        % get minimum and maximum azimuths
        azMax = min(180, angle+abs(tol));
        azMin = max(  0, angle-abs(tol));
        
        % find desired grid points
        id = find( (azHor<=azMax & azHor>=azMin));
        
        % get and sort elevations in transformed coordinate system
        elTrans      = elHor(id);
        [~, idSort]  = sort(elTrans);
        id           = id(idSort);
        
        % get subgrid
        subGrid    = gridIn(id,:);
        
        % transform to inteaural polar
        if exist('sph2hor', 'file')
            [azHor, elHor] = sph2hor(subGrid(:,1), subGrid(:,2));
            subGridHor     = [azHor elHor];
        else
            error('This needs the function ''sph2hor'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
        end
        
    % generate grid
    else
        id         = [];
        
        % get elevation
        el         = (-90:gridIn:270)';
        if el(end) == 270
            el = el(1:end-1);
        end
        
        % construc the grid
        angle   = reshape( angle, numel(angle), 1 );
        subGridHor = [repelem( angle, numel(el), 1 ) repmat( el, numel(angle), 1 )];
        
        % rotate 
        [az, el]  = hor2sph(subGridHor(:,1), subGridHor(:,2));
        [x, y, z] = sph2cart(az/180*pi, el/180*pi, ones(size(az)));
        [az, el]  = cart2sph(y, x, z);
        
        % get spherical coordinates
        subGrid   = round([az el] ./ pi * 180, 4);
        
        % transform to inteaural polar
        if exist('sph2hor', 'file')
            [azHor, elHor] = sph2hor(subGrid(:,1), subGrid(:,2));
            subGridHor     = [azHor elHor];
        else
            error('This needs the function ''sph2hor'' from the SOFA toolbox for Matlab (http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs)')
        end
    end
    
else
    error(['subGrid ''' subGrid '''not existing'])
end

% plot subgrid and input grid
if doPlot
    
    % get only points not included in the subGrid
    if numel(gridIn) > 1
        idIn     = false(numel(gridIn(:,1)), 1);
        idIn(id) = 1;
        [xIn, yIn, zIn] = sph2cart(gridIn(~idIn,1)/180*pi, gridIn(~idIn,2)/180*pi, ones(size(gridIn(~idIn,1))));
    end
    % get points from the subgrid
    [x, y, z] = sph2cart(subGrid(:,1)/180*pi, subGrid(:,2)/180*pi, ones(size(subGrid(:,1))));
    
    % get sphere for nicer plot
    [xS,yS,zS] = sphere(20);
    
    AKf(60,15)
    for n = 1:3
        subplot(1,3,n)
        
        hold on
        if numel(gridIn) > 1
            surf(xS, yS, zS, 'EdgeColor', 'none', 'faceColor', 'w')
            scatter3(xIn, yIn, zIn, '.', 'MarkerFaceColor', [.7 .7 .7], 'MarkerEdgeColor', [.7 .7 .7])
        else
            surf(xS, yS, zS, 'EdgeColor', [.6 .6 .6], 'faceColor', 'w')
        end
        scatter3(x, y, z, '.k')
        
        axis equal; axis off; grid off
        
        if n == 1
            view([90 0])
            title('front view')
        elseif n == 2
            view([-90 0])
            title('back view')
        else
            view([90 90])
            title('top view')
        end
    end
end
