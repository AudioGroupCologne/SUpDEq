% hFigureHandle = AKf(fWidth, fHeight, num)
%
% creates a figure with white background color on center of screen.
% AKf() create a large figure
% AKf(fWidth) creates a square figure
% 
% See AKplotDemo.m for examples
%
% I N P U T:
% fWidth, fHeight - figure widhth and height in cm.
% num (optional)  - number of the figure, If num=fase, the figure will be
%                   invisible
%
% O U T P U T:
% hFigureHandle   - handle to the figure
%
% fabian.brinkmann@tu-berlin.de

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
function hFigureHandle = AKf(fWidth, fHeight, num)

if nargin == 0
    % get size screen size
    if exist('groot', 'file')
        scrn = get(groot, 'ScreenSize');
    else
        scrn = get(0, 'Screensize');
    end
    
    % sometines task bars will shadow the figure - this makes the
    % figure smaller
    try
        % get the size of the matlab main window
        desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
        desktopMainFrame = desktop.getMainFrame;
        desktopDims = desktopMainFrame.getSize;
        % set figure size accordingly
        margin.offset = 10; % works great with Windows 10
        margin.x = scrn(3)-min(scrn(3),desktopDims.getWidth)+margin.offset;
        margin.y = scrn(4)-min(scrn(4),desktopDims.getHeight)+margin.offset;
        margin.X = scrn(3)-2*margin.x;
        margin.Y = scrn(4)-2*margin.y;
        if margin.X > 0 && margin.Y > 0
            scrn = [margin.x margin.y scrn(3)-2*margin.x scrn(4)-2*margin.y];
        else
            scrn = [round(.1*scrn(3)) round(.1*scrn(4)) round(.8*scrn(3)) round(.8*scrn(4))];
        end
    catch
        % if this does not work make the figure smaller by hand
        scrn = [round(.1*scrn(3)) round(.1*scrn(4)) round(.8*scrn(3)) round(.8*scrn(4))];
    end
    
    % create a figure with correct dimensions
    hFigureHandle = figure('outerposition', scrn, 'PaperSize', [scrn(3) scrn(4)], 'PaperPosition', scrn);
    
    % set background color to white
    set(gcf, 'color', [1 1 1]);
    
elseif nargin == 1
    
    hFigureHandle = figure;
    fHeight = fWidth;
    
elseif nargin == 2
    
    hFigureHandle = figure;
    
elseif nargin == 3
    
    if isnumeric(num)
        hFigureHandle = figure(num);
    elseif num
        hFigureHandle = figure;
    else
        hFigureHandle = figure('visible','off');
    end
    
end

if nargin ~= 0
    % get screensize in centimeter
    pixpercm = get(0,'ScreenPixelsPerInch')/2.54;
    screen = get(0,'ScreenSize')/pixpercm;
    % get position for figure
    left   = max([(screen(3)-fWidth)/2 0]);
    bottom = max([(screen(4)-fHeight)/2 0]);
    
    set(hFigureHandle,'PaperUnits', 'centimeters');
    set(hFigureHandle,'Units', 'centimeters');
    
    % paper size for printing
    set(hFigureHandle, 'PaperSize', [fWidth fHeight]);
    % location on printed paper
    set(hFigureHandle,'PaperPosition', [.1 .1 fWidth-.1 fHeight-.1]);
    % location and size on screen
    set(hFigureHandle,'Position', [left bottom fWidth fHeight]);
    
    % set color
    set(hFigureHandle, 'color', [1 1 1])
end

if ~nargout
    clearvars hFigureHandle
end
