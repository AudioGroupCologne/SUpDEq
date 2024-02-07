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
% Visual 3D - sound field data visualisation R13-0306
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% void = sofia_visual3d(mtxData, visualStyle, [colormap])
% ------------------------------------------------------------------------     
% void
% ------------------------------------------------------------------------              
% mtxData      SOFiA 3D-matrix-data [1°steps]
% visualStyle  0 - Globe,      surface colors indicate the intensity
%              1 - Flat,       surface colors indicate the intensity 
%              2 - 3D pattern, extension indicates the intensity
%              3 - 3D pattern, extension+colors indicate the intensity
% color_map    MATLAB colormaps (see MATLAB reference) [optional]
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


function sofia_visual3D(mtxData, visualStyle, color_map)

disp('SOFiA Visual3D Response Visualization R13-0306');

if nargin == 0
   load sofia_exampleMtxData mtxData
end

if nargin < 2
   visualStyle = 0;
end
 
if (size(mtxData,1)~=181 || size(mtxData,2)~=360)
    error('Data type missmatch: SOFiA 3D-matrix-data with 181x360 pixels expected.');
end


set(gcf,'Color', 'w');

mtxData = abs(mtxData);

if visualStyle == 0
    
    if nargin < 3
       color_map = 'Jet'; 
    end
    set(gca, 'NextPlot','add', 'Visible','off');
    axis equal;
    %view(-90,0);
    [x, y, z] = ellipsoid(0, 0, 0, 1, 1, 1);
    globePlot = surf(-x, -y, -z, 'FaceColor', 'none', 'EdgeColor', [1 1 1]);
    alpha = 1;
    set(globePlot, 'FaceColor', 'texturemap', 'CData', mtxData, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    colorbar('location','EastOutside', 'TickLength',[0 0], 'LineWidth', 2, 'FontSize', 16)
    
    view(105,20);  
    
    text(0,0,1.1,'0^oE','fontsize', 14,'color', [.502 .502 .502])
    text(0,0,-1.1,'180^oE','fontsize', 14,'color', [.502 .502 .502])
    text(0,1.1,0,'90^oA','fontsize', 14,'color', [.502 .502 .502])
    text(0,-1.1,0,'270^oA','fontsize', 14,'color', [.502 .502 .502])
    text(-1.1,0,0,'180^oA','fontsize', 14,'color', [.502 .502 .502])
    text(1.1,0,0,'0^oA','fontsize', 14,'color', [.502 .502 .502])
    

elseif visualStyle == 1
    
    if nargin < 3
       color_map = 'Jet'; 
    end
    imagesc(mtxData)
    xlabel('Azimuth in DEG')
    ylabel('Elevation in DEG')

elseif visualStyle == 2
    
    if nargin < 3
       color_map = [0 1 0]; 
    end
    
    mtxData = [mtxData mtxData(:, 1)];
    theta = (0:360)*pi/180;
    phi = (0:180)*pi/180;
    [theta, phi] = meshgrid(theta, phi);
    [X,Z,Y] = sph2cart(theta, phi-pi/2, mtxData);    
    surf(X,Y,Z,'edgecolor','none'); 
    axis equal
    axis off
    box off
    rotate3d on
    light; 
    lighting phong; 
    camzoom(1.3);
       
elseif visualStyle == 3
    
    if nargin < 3
       color_map = 'Jet'; 
    end
    
    mtxData = [mtxData mtxData(:, 1)];
    theta = (0:360)*pi/180;
    phi = (0:180)*pi/180;
    [theta, phi] = meshgrid(theta, phi);
    [X,Z,Y] = sph2cart(theta, phi-pi/2, mtxData);
    set(gca, 'NextPlot','add', 'Visible','off');
    axis equal;
    view(45,45);   
    globePlot = surf(X, Y, -Z, 'FaceColor', 'none', 'EdgeColor', [1 1 1]);
    alpha = 1;
    set(globePlot, 'FaceColor', 'texturemap', 'CData', mtxData, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    colorbar('location','EastOutside', 'TickLength',[0 0], 'LineWidth', 2, 'FontSize', 16)    
    
end

colormap(color_map);