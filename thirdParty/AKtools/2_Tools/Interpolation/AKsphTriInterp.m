% out_vals = AKsphTriInterp(x_azim, x_elev, x, xq_azim, xq_elev, angle_type)
%
%  Triangular interpolation using convex hull triangulation.
%
%  [1] V. Pulkki, Virtual sound source positioning using vector base
%      amplitude panning. J. Audio Eng. Soc. 45(6), 456–466 (1997)
%  [2] F. Zagala, Optimum-phase primal signal and radiation-filter 
%      modelling of musical instruments, Master Thesis (2019)
%
%  Implementation by:   Franz Zotter 2021, zotter@iem.at
%                       Institute of Electronic Music and Acoustics, 
%                       University of Music and Performing Arts Graz
%
%                       David Ackermann, david.ackermann@tu-berlin.de
%                       Audiocommunication Group, TU Berlin
%
%   I N P U T:
%   x_azim, x_elev    - Vectors containing azimuth and elevation of the
%                       sample points. (az=0; el=0) denotes point in front
%                       in front; (90;0) to the left; (0;90) above the
%                       listener
%   x                 - Vector containing the sample points
%   xq_azim,x_q elev  - Vectors containing azimuth and elevation of the
%                       query points
%   angle_type        - 'deg' or 'rad', default: 'deg'
%
%   O U T P U T:
%   out_vals          - interpolated values

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

function out_vals = AKsphTriInterp(x_azim, x_elev, x, xq_azim, xq_elev)

if ~exist('angle_type', 'var')
    angle_type = 'deg';
end

if ~exist('do_plot', 'var')
    do_plot = 0;
end

if(length(x_azim) ~= length(x_elev) || length(x_azim) ~= length(x))
    disp('Error:  Input vectors must have same lengths');
    return;
end

if(strcmpi(angle_type,'deg'))
    x_azim = deg2rad(x_azim);
    x_elev = deg2rad(x_elev);
    xq_azim = deg2rad(xq_azim);
    xq_elev = deg2rad(xq_elev);
end

%% The convex hull of a set of points in 3D space
[X, Y, Z] = sph2cart(x_azim, x_elev, 1);
X = [X,Y,Z];
[Xq, Yq, Zq] = sph2cart(xq_azim, xq_elev, 1);
Xq = [Xq,Yq,Zq];

T = convhull(X);

%% Find assignment of triangles to the 8 octands
oct1=[];
oct2=[];
oct3=[];
oct4=[];
oct5=[];
oct6=[];
oct7=[];
oct8=[];
for tri = 1:size(T,1)
    xpos = (X(T(tri,1),1)>=0)||(X(T(tri,2),1)>=0)||(X(T(tri,3),1)>=0);
    ypos = (X(T(tri,1),2)>=0)||(X(T(tri,2),2)>=0)||(X(T(tri,3),2)>=0);
    zpos = (X(T(tri,1),3)>=0)||(X(T(tri,2),3)>=0)||(X(T(tri,3),3)>=0);
    xneg = (-X(T(tri,1),1)>=0)||(-X(T(tri,2),1)>=0)||(-X(T(tri,3),1)>=0);
    yneg = (-X(T(tri,1),2)>=0)||(-X(T(tri,2),2)>=0)||(-X(T(tri,3),2)>=0);
    zneg = (-X(T(tri,1),3)>=0)||(-X(T(tri,2),3)>=0)||(-X(T(tri,3),3)>=0);
    if  xpos &&  ypos &&  zpos
        oct1=[oct1, tri];
    end
    if  xneg &&  ypos &&  zpos
        oct2=[oct2, tri];
    end
    if  xpos && yneg &&  zpos
        oct3=[oct3, tri];
    end
    if  xpos &&  ypos && zneg
        oct4=[oct4, tri];
    end
    if xneg && yneg &&  zpos
        oct5=[oct5, tri];
    end
    if xneg &&  ypos && zneg
        oct6=[oct6, tri];
    end
    if  xpos && yneg && zneg
        oct7=[oct7, tri];
    end
    if xneg && yneg && zneg
        oct8=[oct8, tri];
    end        
end

%% go through triangles
tol = -1e-3;
xq = nan(size(Xq,1),1);
for q = 1:size(Xq,1)
    xpos = (Xq(q,1)>=0);
    ypos = (Xq(q,2)>=0);
    zpos = (Xq(q,3)>=0);
    if      xpos &&  ypos &&  zpos
        oct=oct1;
    elseif ~xpos &&  ypos &&  zpos
        oct=oct2;
    elseif  xpos && ~ypos &&  zpos
        oct=oct3;
    elseif  xpos &&  ypos && ~zpos
        oct=oct4;
    elseif ~xpos && ~ypos &&  zpos
        oct=oct5;
    elseif ~xpos &&  ypos && ~zpos
        oct=oct6;
    elseif  xpos && ~ypos && ~zpos
        oct=oct7;
    elseif ~xpos && ~ypos && ~zpos
        oct=oct8;
    end
    for tri = oct
        g = (X(T(tri,:),:))' \ Xq(q,:)';
        if sum(g>tol) == 3
            g=g/sum(g);
            xq(q) = x(T(tri,:),:)'*g;
            break
        end
    end
end

out_vals = xq;