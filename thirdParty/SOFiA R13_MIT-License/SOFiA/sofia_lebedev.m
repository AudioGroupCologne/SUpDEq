% /// ASAR/MARA Research Group
%  
% Technology Arts Sciences TH Köln
% Technical University of Berlin
% Deutsche Telekom Laboratories
% University of Rostock
% WDR Westdeutscher Rundfunk
% IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis
% 
% SOFiA Lebedev Grid R13-0306
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% [gridData, Npoints, Nmax] = sofia_lebedev(degree, plot)
% ------------------------------------------------------------------------     
%
% gridData           Lebedev quadrature including weigths(W):
%                    [AZ_1 EL_1 W_1;
%                     AZ_2 EL_2 W_2;
%                     ...
%                     AZ_n EL_n W_n]
%
% Npoints            Total number of nodes
% Nmax               Highest stable grid order  
%
% ------------------------------------------------------------------------
% 
% Order              Lebedev Degree (Number of nodes)
%                    Call sofia_lebedev() to obtain a 
%                    list of valid degrees.
%
% plot               Show a globe plot of the selected grid 
%                    0: Off, 1: On [default]
% 
% This function computes Lebedev quadrature nodes and weigths
% in the SOFiA/VariSphear data format.
% 
%

% CONTACT AND LICENSE INFORMATION:
% 
% /// ASAR/MARA Research Group 
%  
%     [1] Technology Arts Sciences TH Köln
%     [2] Technical University of Berlin 
%     [3] Deutsche Telekom Laboratories 
%     [4] University of Rostock
%     [5] WDR Westdeutscher Rundfunk 
%     [6] IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis toolbox
% 
% Copyright 2011-2017 Benjamin Bernschütz et al.(§)  
% 
% Contact ------------------------------------
% Technology Arts Sciences TH Köln 
% Institute of Communications Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
% 
% phone       +49 221 8275 -2496 
% cell phone  +49 171 4176069 
% mail        rockzentrale 'at' me.com 
% --------------------------------------------
% 
% This file is part of the SOFiA sound field analysis toolbox
%
% Licence Type: MIT License
%
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
%
% (§) Christoph Pörschmann [1]     christoph.poerschmann 'at' th-koeln.de
%     Sascha Spors         [2,3,4] sascha.spors 'at' uni-rostock.de  
%     Stefan Weinzierl     [2]     stefan.weinzierl 'at' tu-berlin.de
%
% External routines for the quadrature calculation in: lebedev_calc.m
% -----------------------------------------------------------------------
% Author: Rob Parrish, The Sherrill Group, CCMST Georgia Tech
% Email : robparrish 'at' gmail.com
% Date  : 03/24/2010
% -----------------------------------------------------------------------
% The external routine is not part of the SOFiA MIT License.


function [gridData, Npoints, Nmax] = sofia_lebedev(degree, plot)

disp('SOFiA Lebedev Grid R13-0306');

degrees_avail=[6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,... 
266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,... 
3074, 3470, 3890, 4334, 4802, 5294, 5810];

if nargin<1
   show_valid_degrees();
   gridData =[];
   Npoints  =[];
   Nmax     =[];
   return
end

if nargin<2
    plot = true;
end

if isempty(find(degrees_avail==degree,1))
   error('Invalid quadrature degree. Call sofia_lebedev() to obtain a list of valid degrees.');
end

% Routines taken from VariSphear waveCapture

[leb] = lebedev_calc(degree);
[theta, phi, notUsed] = cart2sph(leb.x, leb.y, leb.z);
leb.w = leb.w/(4*pi); 
theta = mod(theta, 2*pi);
gridData = [theta, phi+pi/2 ,leb.w];
gridData = sortrows(gridData,2);
gridData = sortrows(gridData,1);

i=1;
turnover=0;
while(1)
    if i>=size(gridData,1)
       break
    end
    c = find(gridData(:,1)==gridData(i,1));
    i = max(c)+1;      
    if turnover == 1 %Toggle
       gridData(c,:)=flipdim(gridData(c,:),1);
       turnover=0;
    else
       turnover=1;
    end
end

Npoints = size(gridData,1);
Nmax    = floor(sqrt(size(gridData,1)/1.3)-1);

if plot
    plot_grid(gridData);
end
end

function show_valid_degrees()
disp(' ');
disp('Valid degrees (number or nodes) are: ');
disp(' ');
disp('6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,')
disp('350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,')
disp('3074, 3470, 3890, 4334, 4802, 5294, 5810')
end

function plot_grid(gridData)

[Xm,Ym,Zm]=sph2cart(gridData(:,1),gridData(:,2)-pi/2,1.01);

colormap Gray;

if size(Xm,1)>1500
    plot3(Xm,Ym,Zm,'marker','.','markerfacecolor','g','color','g','linestyle','none')
else
    plot3(Xm,Ym,Zm,'marker','o','markerfacecolor','g','color','g','linestyle','none')
end
axis off;
hold on;
grid off;
sphere;
axis equal;
rotate3d on;
light;
alpha(.8);
lighting phong;
camzoom(1.4);
hold off;
end

