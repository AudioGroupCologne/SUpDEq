% AKroomSimulationPlot(rs, ISM, sources, ISplot, ISnumber)
% plots the room and reflections paths (optionally) for a room model used
% with AKroomSimulation
%
% Examples:
% AKroomSimulation(rs) plots the room geometry with receiver and source(s)
% AKroomSimulation(rs, ISM) additionally plots all reflection paths
% AKroomSimulation(rs, ISM, 2, 'N', 1) plots first order reflections for
%                                      source two
%
% See AKroomSimulationDemo.m for a use case
% 
% I N P U T
% rs       - strcuct containing the scene description
%            (see AKroomSimulationDemo.m)
% ISM      - output from AKroomSimulation. This is only needed if plotting
%            reflection paths
% sources  - the ID of the sources that should be plotted, eg. [1 3] will
%            plot the first and third source. By default all sources are
%            plotted.
% ISplot   - 'N'  to plot reflections paths of certein image source orders,
%            'id' to plot reflections according to indicees
% ISnumber - scalar or vector that specifies the orders or indices to be
%            plotted
%
% 04/2017  - fabian.brinkmann@tu-berlin.de

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
function AKroomSimulationPlot(rs, ISM, sources, ISplot, ISnumber)

%% ------------------------------------------------------ default parameter
if ~exist('sources', 'var')
    sources = 1:size(rs.srcPos, 1);
end

if ~exist('ISM', 'var')
    ISM = [];
end

if ~isempty(ISM)
    if ~exist('ISplot', 'var')
        ISplot = N;
    end
    if ~exist('ISnumber', 'var')
        ISnumber = 0;
    end
end

if exist('ISplot', 'var')
    if strcmpi(ISplot, 'N')
        ISnumber = min(ISnumber, 3);
    end
end

% --- check if the source view was passed --- %
if isempty(rs.srcView)
    
    % get source position in spherical coordinates
    if strcmpi(rs.srcCoordinates, 'spherical')
        tmp = rs.srcPos;
    else
        [az, el] = cart2sph(rs.srcPos(:,1)-rs.recPos(1), rs.srcPos(:,2)-rs.recPos(2), rs.srcPos(:,3)-rs.recPos(3));
        tmp      = [az el];
        tmp      = tmp*180/pi;
    end
    
    % wrap to 360 deg
    tmp(:,1) = mod(tmp(:,1), 360);
    
    % mirror to get source rotation
    rs.srcView = [mod(180+tmp(:,1), 360) -tmp(:,2)];
    clear az el tmp
    
end

% --- get the source position(s) in carthesian coordinates --- %
if strcmpi(rs.srcCoordinates, 'spherical')
    [x, y, z] = sph2cart(rs.srcPos(:,1)/180*pi, rs.srcPos(:,2)/180*pi, rs.srcPos(:,3));
    rs.srcPos = [x + rs.recPos(1) y + rs.recPos(2) z + rs.recPos(3)];
    
    rs.srcCoordinates = 'carthesian';
    
    clear x y z
end

%% --------------------------------------------------------------- plotting
AKf(30,20)
hold on

% -------- plot the room -------- %
patch([0 0 0 0], [0 rs.L(2) rs.L(2) 0], [0 0 rs.L(3) rs.L(3)], 'k', 'FaceAlpha', .1)
patch([rs.L(1) rs.L(1) rs.L(1) rs.L(1)], [0 rs.L(2) rs.L(2) 0], [0 0 rs.L(3) rs.L(3)], 'k', 'FaceAlpha', .1)
patch([0 rs.L(1) rs.L(1) 0], [0 0 0 0], [0 0 rs.L(3) rs.L(3)], 'k', 'FaceAlpha', .1)
patch([0 rs.L(1) rs.L(1) 0], [rs.L(2) rs.L(2) rs.L(2) rs.L(2)], [0 0 rs.L(3) rs.L(3)], 'k', 'FaceAlpha', .1)
patch([0 rs.L(1) rs.L(1) 0], [0 0 rs.L(2) rs.L(2)], [0 0 0 0], 'k', 'FaceAlpha', .1)
patch([0 rs.L(1) rs.L(1) 0], [0 0 rs.L(2) rs.L(2)], [rs.L(3) rs.L(3) rs.L(3) rs.L(3)], 'k', 'FaceAlpha', .1)


% -------- plot the rays -------- %
if ~isempty(ISM)
    
    for mm = 1:numel(sources)
        
        % source and receiver position
        R = rs.recPos;
        S = rs.srcPos(sources(mm),:);
        
        % get the ray information
        wallLog = ISM(sources(mm)).wallLog;
        Nis     = zeros(size(wallLog,1),1);
        for nn = 1:size(wallLog,1)
            Nis(nn) = ISM(sources(mm)).wallLog(nn).N;
        end
        
        % search rays of desired order
        if strcmpi(ISplot, 'N')
            tmp = ISnumber;
            IS_id = find(Nis == tmp(1));
            for nn = 2:numel(tmp)
                IS_id = [IS_id; find(Nis == tmp(nn))]; %#ok<AGROW>
            end
        end
        
        colors = 'pgy';
        
        % loop accros rays to be plotted
        for nn = 1:numel(IS_id)
            
            % current color
            cc = mod( Nis(IS_id(nn)), 6 ) + 1;
            cc = AKcolors(colors(cc));
            
            % plot the ray
            if isempty(wallLog(IS_id(nn)).ID)
                plot3([S(1) R(1)], [S(2) R(2)], [S(3) R(3)], 'color', cc, 'lineWidth', 2)
            else
                plot3([R(1) wallLog(IS_id(nn)).position(1,1)], [R(2) wallLog(IS_id(nn)).position(1,2)], [R(3) wallLog(IS_id(nn)).position(1,3)], 'color', cc, 'lineWidth', 2)
                
                for pp = 2:numel(wallLog(IS_id(nn)).ID)
                    plot3([wallLog(IS_id(nn)).position(pp-1,1) wallLog(IS_id(nn)).position(pp,1)], ...
                        [wallLog(IS_id(nn)).position(pp-1,2) wallLog(IS_id(nn)).position(pp,2)], ...
                        [wallLog(IS_id(nn)).position(pp-1,3) wallLog(IS_id(nn)).position(pp,3)], 'color', cc, 'lineWidth', 2)
                end
                
                plot3([S(1) wallLog(IS_id(nn)).position(end,1)], [S(2) wallLog(IS_id(nn)).position(end,2)], [S(3) wallLog(IS_id(nn)).position(end,3)], 'color', cc, 'lineWidth', 2)
            end
        end
    end
end

% -------- plot receiver and source(s) -------- %
% get a sphere for marking the positions
r            = min(1, .05 * min(rs.L));
[xS, yS, zS] = sphere(10);
xS = xS * r;
yS = yS * r;
zS = zS * r;

% plot the receiver position
surf(xS+rs.recPos(1), yS+rs.recPos(2), zS+rs.recPos(3), 'FaceColor', AKcolors('r'), 'EdgeColor', 'none')

% plot the receiver orientation
recAz = rs.recView(1);
recEl = rs.recView(2);
[x, y, z] = sph2cart(mean(recAz)/180*pi, mean(recEl)/180*pi, 1.2*r);
plot3([rs.recPos(1) x+rs.recPos(1)], [rs.recPos(2) y+rs.recPos(2)], [rs.recPos(3) z+rs.recPos(3)], 'k', 'LineWidth', 2)

% plot the source position(s)
for nn = 1:numel(sources)
    
    x = rs.srcPos(sources(nn),1);
    y = rs.srcPos(sources(nn),2);
    z = rs.srcPos(sources(nn),3);
    
    % plot the source
    surf(xS+x, yS+y, zS+z, 'FaceColor', AKcolors('b'), 'EdgeColor', 'none')
    
    % plot the source orientation
    srcAz = rs.srcView(sources(nn),1);
    srcEl = rs.srcView(sources(nn),2);
    [xR, yR, zR] = sph2cart(srcAz/180*pi, srcEl/180*pi, 1.2*r);
    plot3([x x+xR], [y y+yR], [z z+zR], 'k', 'LineWidth', 2)
    
end

% -------- format the plot -------- %
view([-45 15])
axis equal
xlabel x; ylabel y; zlabel z
set(gca, 'xTick', [], 'yTick', [], 'zTick', [])
rotate3d on
title({'Scene geometry (source = blue, receiver = red)' 'Black dots indicate viewing direction'})