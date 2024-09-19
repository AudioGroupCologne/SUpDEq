% out_vals = AKsphSplineInterp(x_azim, x_elev, x, xq_azim, xq_elev, m, lambda, angle_type, do_plot)
%
%  Spherical interpolation and smoothing using thin plate splines.
%  See: AKsphSplineInterpDemo.m for examples
%
%  [1] Grace Wahba: "Spline Interpolation and Smoothing on the Sphere."
%      SIAM J. Sci. Stat. Comput., 2(1):5-16, 1981
%  [2] Grace Wahba: "Erratum: Spline Interpolation and Smoothing on the
%      Sphere." SIAM J. Sci. Stat. Comput., 2(2):385-386, 1982
%
%  Implementation by: Tobias Muenzer 2014, txmunzerx@gmail.com
%  Audiocommunication Group, TU Berlin
%
%   I N P U T:
%   x_azim, x_elev    - Vectors containing azimuth and elevation of the
%                       sample points. (az=0; el=0) denotes point in front
%                       in front; (90;0) to the left; (0;90) above the
%                       listener
%   x                 - Vector containing the sample points
%   xq_azim,x_q elev  - Vectors containing azimuth and elevation of the
%                       query points
%   m                 - Order of the interpolation m = {1, 2, 3}, default: 1
%   lambda            - Smoothing factor. Higher values result in stronger
%                       smoothing default: 0
%   angle_type        - 'deg' or 'rad', default: 'deg'
%   do_plot           - Plot interpolation results. Default: 0
%                       0: Do not plot
%                       1: On azimuth-elevation plane
%                       2: 3d plot on sphere
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
function out_vals = AKsphSplineInterp(x_azim,x_elev,x, xq_azim, xq_elev,m,lambda,angle_type,do_plot)

if ~exist('angle_type', 'var')
    angle_type = 'deg';
end
if ~exist('m', 'var')
    m = 1;
end
if ~exist('lambda', 'var')
    lambda = 0;
end
if ~exist('do_plot', 'var')
    do_plot = 0;
end

if(length(x_azim) ~= length(x_elev) || length(x_azim) ~= length(x))
    disp('Error:  Input vectors must have same lengths');
    return;
end

if(length(xq_azim) ~= length(xq_elev))
    disp('Error:  Query vectors must have same lengths');
    return;
end


if(strcmpi(angle_type,'deg'))
    x_azim = deg2rad(x_azim);
    x_elev = deg2rad(x_elev);
    xq_azim = deg2rad(xq_azim);
    xq_elev = deg2rad(xq_elev);
end

%% Solve the equation (R_n + n*lambda*I)*c + dT = z  T'c = 0 for c and d
%  where z are the values at the sample points and R_n is a matrix with
%  the m,n th entry defined by R(P_m,P_n)

n = length(x_azim);

%% A Matrix
A = zeros(n+1,n+1);
for i=1:n
    A(1:n,i) = R_function(x_azim,x_elev,x_azim(i),x_elev(i),m);
end
A(:,end) = ones(n+1,1);
A(end,:) = ones(1,n+1);

%% smoothing factor
I = zeros(n+1,n+1);
I(1:end-1,1:end-1) = eye(n);
A = A + n*lambda*I;

%% B Vector
B = [x;0];

%% Solve A*c = B
c = linsolve(A,B);


%% Interpolate points
%  Calculate Sum i=1:n(c_i(R(P,P_i)) + d

u = 0;
for i=1:n
    u = u + c(i)*R_function(x_azim(i),x_elev(i),xq_azim,xq_elev,m);
end
u = u + c(n+1); % add d

out_vals = u;

if(do_plot == 1)
    plotPoints2d(x_azim,x_elev,x, xq_azim, xq_elev, out_vals);
end
if(do_plot == 2)
    plotPoints3d(x_azim,x_elev, xq_azim, xq_elev, out_vals);
end

end


function q = R_function(phi1, theta1, phi2,theta2,m)

% %angle between sampel points and query points
% g = acos(sin(theta1).*sin(theta2)  + cos(theta1).*cos(theta2).*cos(phi1-phi2));
% 
% % defined in [1],eq. 3.3
% z = cos(g);
% 
% %Equations for m = 1..3 defined in [1],table 1 (haversine function with z = cos(g))
% w = (1 - z)./2;

% if the two points are close together one might end up with cos(x) = 0.99999999, 
% leading to an inaccurate answer. Since the haversine formula uses sines, 
% it avoids that problem.
w = sin(0.5*(theta1-theta2)).^2 + cos(theta1).*cos(theta2).*sin(0.5*(phi1-phi2)).^2;

a = log(1+(1./sqrt(w)));
c = 2*sqrt(w);

lim_val = 0;    %value of polynom for d->0

if(m == 1)
    k = 3/2;
    q = 2* a.*w - c+1;
    q = (q/factorial(2*k-2) - 1/factorial(2*k-1))./(2*pi);
    
    lim_val = 1;
    lim_val = (lim_val/factorial(2*k-2) - 1/factorial(2*k-1))./(2*pi);
elseif(m == 2)
    k = 2;
    q = (a.*(12.*w.^2-4.*w)-6.*c.*w+6.*w+1)/2;
    q = (q/factorial(2*k-2) - 1/factorial(2*k-1))./(2*pi);

    lim_val = 1/2;
    lim_val = (lim_val/factorial(2*k-2) - 1/factorial(2*k-1))./(2*pi);
elseif(m == 3)
    k = 5/2;
    q = (a.*(60*w.^3-36*w.^2)+30*w.^2+c.*(8*w-30*w.^2)-3*w+1)/3;
    q = (q/factorial(2*k-2) - 1/factorial(2*k-1))./(2*pi);
    
    lim_val = 1/3;
    lim_val = (lim_val/factorial(2*k-2) - 1/factorial(2*k-1))./(2*pi);
end

f = isnan(q);
q(f) = lim_val;

end


function plotPoints2d(x_azim,x_elev,x, xq_azim, xq_elev, xq)

%% PLOT ON 2D SURFACE

n = length(x);
nq = length(xq);

pointSize = [6*ones(1,nq) 6*ones(1,n)];

all_azim = [xq_azim;x_azim] /pi*180;
all_elev = [xq_elev; x_elev]/pi*180;
all_data = [xq;x];

if isempty(findall(0,'Type','Figure'))
    AKf(20,10)
end
scatter(all_azim,all_elev,pointSize,all_data,'filled', 'SizeData', 30);
hold on;
plot(x_azim/pi*180,x_elev/pi*180,'xk');
hold off;
view([0 90])
xlabel('Azimuth [deg]');
ylabel('Elevation [deg]');
title('Spherical interpolation results');
axis([0 360 -90 90]);

end

function plotPoints3d(x_azim,x_elev, xq_azim, xq_elev, xq)

%% PLOT ON 3D SURFACE

[X, Y, Z] =    sph2cart(x_azim,x_elev,1);
[Xq, Yq, Zq] = sph2cart(xq_azim,xq_elev,1);

if isempty(findall(0,'Type','Figure'))
    AKf(20,10)
end
scatter3(Xq,Yq,Zq,6,xq,'filled', 'SizeData', 30);
hold on;
plot3(X,Y,Z,'xk');
hold off;
axis equal

title('Spherical interpolation results');

end
