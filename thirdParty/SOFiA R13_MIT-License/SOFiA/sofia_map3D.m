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
% Map3D - Sound Field Visualisation based on Map Projections R12-0104
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% void = sofia_map3d(mtxData, projection, textSize)
% ------------------------------------------------------------------------
% void         No return.
% ------------------------------------------------------------------------
% mtxData      SOFiA 3D-matrix-data [1°steps]
% projection   Map Projection Type: 
%               'mollweid': Mollweide Projection (#default) 
%               'hammer'  : Hammer Projection
%               'aitoff'  : Hammer/Aitoff Projection
%
% Important: MATLAB Mapping Toolbox Required 
%            http://www.mathworks.de/products/mapping/
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

function sofia_map3D(mtxData, projection)

disp('SOFiA Map3D Visual Map Projection R13-0306');

    mappingToolboxExists = false;
    toolboxes = ver;
    for i = 1: size(toolboxes,2)
         if ~isempty(strfind(toolboxes(i).Name,'Mapping Toolbox'))
             disp('>>> Mapping toolbox: ok');
             mappingToolboxExists = true;
             break
         end
    end
if ~mappingToolboxExists
    disp('>>> Mapping toolbox: failed');
    error('The MATLAB Mapping Toolbox is required and could not be found.');
end

if nargin == 0
   load sofia_exampleMtxData mtxData
end

if nargin < 2
    projection = 'mollweid';
end

if ~strcmp(projection,'mollweid') && ~strcmp(projection,'hammer') && ~strcmp(projection,'aitoff')
error('Type of map projection is not valid');
end

if (size(mtxData,1)~=181 || size(mtxData,2)~=360)
    error('Data type missmatch: SOFiA 3D-matrix-data with 181x360 pixels expected.');
end

mtxData = abs(mtxData);

clf();
set(gcf,'Color', 'w');
set(gcf, 'units', 'centimeters', 'pos', [2 3 20 10])

mtxData = flipud(mtxData);
mtxData = fliplr(mtxData);
reverse = 1;   % Vertical Label Location Left(0)/Right(1)
textSize = 12; % Size of the Map Labels


% Create referencing matrix
mtxDataR = makerefmat('RasterSize', size(mtxData), ...
    'Latlim', [-90 90], 'Lonlim', [0 360]);

% Set up proj
mapaxes = axesm('MapProjection',projection,'Frame','off','Grid','on','ParallelLabel','on', 'MeridianLabel','on');
if reverse == 0
    setm(gca, 'PLabelMeridian', -90, 'PlineLocation', 30)
else
    setm(gca, 'PLabelMeridian', 90, 'PlineLocation', 30)
end

setm(gca, 'MLabelParallel', 0, 'MlineLocation', 30)

setm(gca, 'FontColor', [0.7, 0.7, 0.7], 'FontSize', textSize, 'FontWeight', 'bold')
setm(gca, 'LabelFormat','none', 'LabelRotation','off', 'maplonlimit', [-180 180]);
spacing = [200 400];
setm(mapaxes, 'glinewidth', 1, 'glinestyle', ':', 'gcolor', [0.7, 0.7, 0.7])

% Display data mapped to the graticule
meshm(mtxData, mtxDataR, spacing);
h=findobj(gcf, 'type', 'text');

if strcmp(projection,'hammer') || strcmp(projection,'aitoff') 
    
    if textSize <= 30
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
    set(h(3),'String','30^{\circ}');
    set(h(4),'String','');
    set(h(5),'String','60^{\circ}');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','120^{\circ}');
    set(h(10),'String','');
    set(h(11),'String','150^{\circ}');
    set(h(12),'String','');
    set(h(13),'String','');
    
    %Parallel
    set(h(14),'String','210^{\circ}');
    set(h(15),'String','240^{\circ}');
    set(h(16),'String','270^{\circ}');
    set(h(17),'String','300^{\circ}');
    set(h(18),'String','330^{\circ}');
    
    set(h(20),'String','30^{\circ}');
    set(h(21),'String','60^{\circ}');
    set(h(22),'String','90^{\circ}');
    set(h(23),'String','120^{\circ}');
    set(h(24),'String','150^{\circ}');
    
    else
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
   % set(h(3),'String','30^{\circ}');
    set(h(3),'String','');
    set(h(4),'String','');
    set(h(5),'String','');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','');
    set(h(10),'String','');
   % set(h(11),'String','150^{\circ}');
    set(h(11),'String','');
    set(h(12),'String','');
    set(h(13),'String','');
    
    if reverse == 0
        set(h(16),'String','270^{\circ}');
        set(h(22),'String','90^{\circ}');
    else
        %set(h(22),'String','270^{\circ}');
        %set(h(16),'String','90^{\circ}');
        set(h(22),'String','');
        set(h(16),'String','');
        set(h(19),'String','');    
    end
    
    %Parallel
    set(h(14),'String','');
    set(h(15),'String','');    
    set(h(17),'String','');
    set(h(18),'String','');    
    set(h(20),'String','');
    set(h(21),'String','');   
    set(h(23),'String','');
    set(h(24),'String','');
    end
        
elseif strcmp(projection,'mollweid') 
   
    if textSize <= 30   
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
    set(h(3),'String','30^{\circ}');
    set(h(4),'String','');
    set(h(5),'String','60^{\circ}');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','120^{\circ}');
    set(h(10),'String','');
    set(h(11),'String','150^{\circ}');
    set(h(12),'String','');
    set(h(13),'String','');
    
    %Parallel
    set(h(14),'String','');
    set(h(15),'String','210^{\circ}');
    set(h(16),'String','240^{\circ}');
    set(h(17),'String','270^{\circ}');
    set(h(18),'String','300^{\circ}');
    set(h(19),'String','330^{\circ}');
    
    set(h(20),'String','0^{\circ}');
    set(h(21),'String','30^{\circ}');
    set(h(22),'String','60^{\circ}');
    set(h(23),'String','90^{\circ}');
    set(h(24),'String','120^{\circ}');
    set(h(25),'String','150^{\circ}');
    set(h(26),'String','');
    
    
    else
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
    set(h(3),'String','30^{\circ}');
    set(h(4),'String','');
    set(h(5),'String','');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','');
    set(h(10),'String','');
    set(h(11),'String','150^{\circ}');
    set(h(12),'String','');
    set(h(13),'String','');
    
    %Parallel
    set(h(14),'String','');
    set(h(15),'String','');
    set(h(16),'String','');
    if reverse == 0
        set(h(17),'String','270^{\circ}');
        set(h(23),'String','90^{\circ}');
    else
        set(h(23),'String','270^{\circ}');
        set(h(17),'String','90^{\circ}');
    end
    set(h(20),'String','0^{\circ}');
    set(h(18),'String','');
    set(h(19),'String','');   
    set(h(21),'String','');
    set(h(22),'String','');    
    set(h(24),'String','');
    set(h(25),'String','');
    set(h(26),'String','');        
    end
        
end

axis off
box off
colorbar();



