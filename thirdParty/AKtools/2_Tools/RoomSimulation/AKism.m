function [d, A, Raz, Rel, Saz, Sel, X, Xrel, N, wallLog] = AKism(L, S, Sview, R, Rview, alpha, truncation, c)
% Calculates image source model (ISM) for a cubic rooms with real, and
% frequency dependend reflection coefficients for each wall.
%
% For example usage see
% AKroomSimulationDemo.m and AKroomSimulation.m
%
% I N P U T
% L          - array containing room dimensions in m [x y z]. One corner of
%              the room is set in the origin of co-ordinates. The room lies
%              in positive x,y,z-direction.
% S          - three element vector containing source position in m [x y z]
%              (all positive)
% Sview      - two element vector specifying the source view in
%              [azimuth elevation]
%              azimuth: 0 deg. = front, 90 deg. = left, etc.
%              (front defined by positive x-axis, left by positive y-axis)
%              elevation: 0 deg. = front, 90 deg. = top, -90 deg. = bottom.
%              (top defined by positive z-axis)
% R          - array containing receiver position in m [x y z], all
%              positive
% Rview      - two element vector specifying the receiver view in
%              [azimuth elevation] (see Sview for coordinate convention)
% alpha      - absorption coefficients:
%              [x1,1 x2,1 y1,1 y2,1 z1,1 z2,1
%               x1,2 x2,2 y1,2 y2,2 z1,2 z2,2]
%                 .     .    .    .    .    .
%                 .     .    .    .    .    .
%                 .     .    .    .    .    .
%               x1,K x2,K y1,K y2,K z1,K z2,K],
%              where x, y, z denote walls with constant x, y, z coordinates
%              , the first index denotes the wall position (1 beeing closer
%              to the origin of coordinates) and the second index the
%              frequency
% truncation - specifies the truncation of the image source model
%                : to include reflection up to order M pass
%                  {'N' M}
%                : to include all refluction up to T seconds pass
%                  {'t' T}
%                : to truncate after P times the estimated perceptual
%                  mixing time pass:
%                  {'tmix' P}
% c          - speed of sound in m/s (default = 343)
%
%
% O U T P U T
% d       - distance of image sources to the receiver
% A       - amplitude of image sources
% Raz     - azimuth incident angle at the receiver
%           (see Sview for coordinate convention)
% Rel     - elevation incident angle at the receiver
%           (see Sview for coordinate convention)
% Saz     - azimuth exit angle at the source
% Sel     - elevation exit angle at the source
% X       - absoulte image source coordinates
% Xrel    - image source coordinates relative to the receiver
% N       - image source order (i.e. number of refelctions)
% wallLog - struct with M entries containing the fields 'ID' and 'position'
%         'ID'        gives the wall(s) that were hit in reverse order,i.e.
%                     the wall that was hit last is listed first (ID are
%                     numbered in the order specified by alpha, i.e. ID1 
%                     reffers to wall x1, ID2 to x2, ID3 to y1, etc.)
%         'position'  gives the intersection points with the walls in
%                     x/y/z coordinates
%         'N'         gives the image source order
%         NOTE that walls is only calculated up to the third image source
%         order and is intended for visualzation only.
%
%
% [1] Lehmann, E. A. and Johansson, A. M. (2008): "Prediciton of energy
%     decay in room impulse responses simulated with an image-source model."
%     In: J. Acoust. Soc. Am. 124(1):269-277.
% [2] Brinkmann, F. & Erbes, V. & Weinzierl, S.: "Extending the closed form
%     image source model for source directivty." Fortschritte der Akustik -
%     DAGA 2018, Munich, Germany, March 2018.
%
% fabian.brinkmann@tu-berlin.de
% 11/2015   - initial code
% 04/2017   - added truncation criteria
% 10/2017   - added the calculation of the source exit angles

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

% set default parameter
if ~exist('c', 'var')
    c = 343;
end

% process the truncation parameter, and set limits for the loops
if ~any(strcmpi(truncation{1}, {'N' 't' 'tmix'}))
    error('AKism:input', ['''' truncation{1} ' is an invalid value for the truncation parameter'])
end

if strcmpi(truncation{1}, 'N')
    if truncation{2} == 0
        Nx = 0;
        Ny = 0;
        Nz = 0;
    else
        Nx = 1-truncation{2}:truncation{2};
        Ny = Nx;
        Nz = Nx;
    end
    
    % maximum allowed IS distance (not needed if truncation is according to N)
    dMax = inf;
    
end
    
if strcmpi(truncation{1}, 'tmix')
    tmix          = 20 * prod(L)/(2*L(2)*L(3) + 2*L(1)*L(3) + 2*L(1)*L(2)) + 12;
    truncation{2} = truncation{2} * tmix/1000;
    
    clear tmix
end
    
if any( strcmpi(truncation{1}, {'t' 'tmix'}) )
    
    % convert maximum IS delay to maximum IS distance
    dMax = truncation{2}*c;
    
    % add the source-to-receiver distance
    % (mixing time starts at the direct sound)
    if strcmpi(truncation{1}, 'tmix')
        dMax = dMax + sqrt( sum( (S-R).^2 ) );
    end
    
    Nx = ceil( dMax / ( 2*L(1) ) );
    Nx = -Nx:Nx;
    
    Ny = ceil( dMax / ( 2*L(2) ) );
    Ny = -Ny:Ny;
    
    Nz = ceil( dMax / ( 2*L(3) ) );
    Nz = -Nz:Nz;
end

% get refelction coefficients
b = sqrt(1-alpha);

% number of image sources
M = numel(Nx) * numel(Ny) * numel(Nz) * 2^3;

% allocate memory
Xrel = zeros(M, 3);         % relatice IS position
X    = Xrel;                % absolute IS position
d    = zeros(M, 1);         % distance between IS and receiver
A    = zeros(M, size(b,1)); % IS amplitude
hits = zeros(M, 6);         % number of hits per wall
N    = zeros(M, 1);         % image source order

% ------------------------------------------------------- get image sources
p = 1;  % counter for image sources
for u = 0:1
    for v = 0:1
        for w = 0:1
            for l = Nx
                for m = Ny
                    for n = Nz
                        
                        % wall hits (not chronologically ordered!)
                        hits(p,:) = [abs(u-l) abs(l) abs(v-m) abs(m) abs(w-n) abs(n)];
                        
                        % image source order
                        N(p) = sum(hits(p,:), 2);
                        
                        % check if we need to calculate the remaining parameters
                        if ~strcmpi(truncation{1}, 'N') || ( strcmpi(truncation{1}, 'N') && N(p) <= truncation{2} )
                            
                            % Absolute position of current (image) source
                            X(p,:)    = [(1-2*u)*S(1)+2*l*L(1) ...
                                         (1-2*v)*S(2)+2*m*L(2) ...
                                         (1-2*w)*S(3)+2*n*L(3)];
                            
                            % position of (image) source relative to the receiver
                            Xrel(p,:) = X(p,:) - R;
                            
                            % distance to receiver (cf. [1], eq. 8)
                            d(p) = norm(Xrel(p,:));
                            
                            % check if we need to calculate the remaining parameters
                            if d(p) <= dMax
                                
                                % amplitude calculated from reflection coefficients b
                                % (cf. [1], eq. 6)
                                A(p,:) = (b(:,1).^abs(u-l) .* b(:,2).^abs(l) .*...
                                          b(:,3).^abs(v-m) .* b(:,4).^abs(m) .*...
                                          b(:,5).^abs(w-n) .* b(:,6).^abs(n)) / ...
                                          (4*pi*d(p));
                                
                                % increment counter
                                p = p+1;
                            end
                            
                        end
                        
                    end
                end
            end
        end
    end
end

clearvars u v w l n m Nx Ny Nz

% ------------------------------------------------------------- sort output

M = p-1;

[d, id] = sort(d(1:M));
A       = A(id,:);
X       = X(id,:);
Xrel    = Xrel(id,:);
hits    = hits(id,:);
N       = N(id);

clear id p

% ----------------------------------- incident sound angles at the receiver

% angles for neutral receiver view
Raz = atan2d(Xrel(:,2), Xrel(:,1));
Rel = 90 - acosd(Xrel(:,3)./d);

% rotate according to receiver view
if any(Rview)
    [Raz, Rel] = AKroomSimulationRotation(Raz, Rel, Rview(1), Rview(2));
end

Raz = mod(Raz, 360);

% ----------------------------------------- exit sound angles at the source

% --- angles for neutral source view
% vector pointing from receiver to image sources
Xp = -Xrel;

% number of wall hits on walls with constant x,y, and z-coordinates
id = [sum(hits(:,1:2), 2) sum(hits(:,3:4), 2) sum(hits(:,5:6), 2)];
% indicate if number of walls hits are odd
id = mod(id, 2);
id = logical( id );
% mirror the source vector in case of odd number of wall hits
Xp(id) = -Xp(id);

% get the exit angles
Saz = atan2d( Xp(:,2), Xp(:,1) );
Sel = 90 - acosd( (Xp(:,3)) ./ d );

% --- rotate according to source view
if any(Sview)
    [Saz, Sel] = AKroomSimulationRotation(Saz, Sel, Sview(1), Sview(2));
end


% --------------------------------------------------- backwards ray tracing
% I N P U T
% X    - [M x 3] array containing absolute position of source and image
%        source(s) in carthesian coordinates [x1 y1 z1; x2 y2 z2; ...]
% R    - Reiceiver position in cartesian coordinates [x y z]
% L    - Room dimension in x/y/z direction [x y z]
% hits - [M x 6] matrix that holds the number of reflections on each wall
%         wall ID    1    2    3    4    5    6
%                 [x1,1 x2,1 y1,1 y2,1 z1,1 z2,1
%                  x1,2 x2,2 y1,2 y2,2 z1,2 z2,2]
%                    .     .    .    .    .    .
%                    .     .    .    .    .    .
%                    .     .    .    .    .    .
%                  x1,M x2,M y1,M y2,M z1,M z2,M],
%        where x, y, z denote walls with constant x, y, z coordinates, the
%        first index denotes the wall position (1 beeing closer to the
%        origin of coordinates) and the second index the frequency
%
%
% O U T P U T
% wallLog - struct with M entries containing the fields 'ID' and 'position'
%           'ID'        gives the wall(s) that were hit (see input 'hits')
%                       in reverse order (i.e. the wall that was hit last
%                       is listed first)
%           'position'  gives the intersection points with the walls in
%                       x/y/z coordinates

% generate struct to save the walls that were hit and the position, where
% the walls were hit
idN3    = find(N<=3);
wallLog = struct('ID', [], 'position', [], 'N', 0);
wallLog = repmat(wallLog, [numel(idN3) 1]);

% allocate memory to save the first point of reflection to calculate the
% source exit angle
Xref = zeros(size(X));
% for the direct sound we use the position of the receiver
Xref(1,:) = R;

% normal wall vectors (for walls 1 to 6 ordered according to alpha)
wNormal = [ 1  0  0
           -1  0  0
            0  1  0
            0 -1  0
            0  0  1
            0  0 -1];

% an arbitrary point on the wall (for walls 1 to 6 ordered according to alpha)
wPoint = [0    0    0
          L(1) 0    0
          0    0    0
          0    L(2) 0
          0    0    0
          0    0    L(3)];
      
% coordinate (1:x, 2:y, 3:z) and value for mirroring (for walls 1 to 6 ordered according to alpha)
wMirror = [1 0
           1 L(1)
           2 0
           2 L(2)
           3 0
           3 L(3)];

for ll = 1:numel(idN3)
    
    nn = idN3(ll);
    
    % allocate space
    wallLog(ll).ID       = zeros(N(nn), 1);
    wallLog(ll).position = zeros(N(nn), 3);
    wallLog(ll).N        = N(nn);
    
    % starting and end position of the ray
    posStart = X(nn,:);
    posEnd   = R;
    
    % memorize the latest wall three wall hits to make sure that rays don't
    % get caught in edges or corners
    lastWall = [0 0 0];
        
    % trace back the path
    for mm = 1:N(nn)
        
        % possible walls for reflections
        walls = find(hits(nn,:) ~= 0);
        % same wall can't be hit two times in a row
        walls = walls(walls ~= lastWall(1));
        % check if we hit an edge or corner:
        % latest two or three intersection points would be identical in
        % this case
        if mm > 3
            dI = diff( wallLog(ll).position(mm-3:mm-1,:) );
            % check for corner
            if all( abs(dI) < 2e-6 )
                walls = walls(walls ~= lastWall(2) & walls ~= lastWall(3));
            % check for edge
            elseif all( abs(dI(2,:)) < 2e-6 )
                walls = walls(walls ~= lastWall(2));
            end
        elseif mm > 2
            dI = diff( wallLog(ll).position(mm-2:mm-1,:) );
             % check for edge
            if all( abs(dI) < 2e-6 )
                walls = walls(walls ~= lastWall(2));
            end
        end
        
        for ww = 1:numel(walls)
            
            currWall = walls(ww);
            
            % check if the intersection with the wall exists
            intersect = AKwallIntersect( wNormal(currWall,:), wPoint(currWall,:), posStart, posEnd, L );
            
            % save intersection point if it exists
            if ~islogical(intersect)
                
                % save ID of the wall that was hit and the intersection
                wallLog(ll).ID(mm)         = currWall;
                wallLog(ll).position(mm,:) = intersect;
                
                if mm == N(nn)
                    Xref(nn,:) = intersect;
                end
                
                % distance of current IS to current wall
                currDistW = wMirror(currWall, 2) - posStart(wMirror(currWall));
                % get the start point by mirroring on the current wall
                posStart(wMirror(currWall, 1)) = posStart(wMirror(currWall, 1)) + 2*currDistW;
                % get the new end point
                posEnd = intersect;
                
                % decrease hit list
                hits(nn, currWall) = hits(nn, currWall) - 1;
                
                % memoraize latest wall hit
                lastWall(2:3) = lastWall(1:2);
                lastWall(1)   = currWall;
                
                break
                
            end
            
        end
    end
end