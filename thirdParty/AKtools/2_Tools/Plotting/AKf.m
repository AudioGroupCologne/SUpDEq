% hFigureHandle = AKf(fWidth, fHeight, num)
%
% creates a figure with white background color on center of screen.
% AKf() create a maximized figure
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
% 
% v1         Fabian Brinkmann, fabian.brinkmann@tu-berlin.de,
%            Audio Communicatin Group, TU Berlin,
%            initial dev
% v2         Hannes Helmholz, helmholz@campus.tu-berlin.de,
%            Audio Communicatin Group, TU Berlin,
%            introduction of AKf() creating a figure with somehow maximum
%            screen dimensions
% v2 04/2018 modification of AKf() to use MAXIMIZE for creating a properly
%            maximized figure window

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

    hFigureHandle = figure;

    % set background color to white
    set(gcf, 'color', [1 1 1]);
    
    % The MAXIMIZE function throws a warning about the used API being
    % depricated in the future. The warning occured but functions seem to
    % work fine: R2017b (9.3).
    %
    % Other versions have not been tested yet. Since what version does it
    % occur? From what version will the API be fully discharged? The
    % warning surpression with verLessThan() should be adjusted or
    %  new relaeses of maximize() should be introduced!
    %
    % To not clog up the command windows this warning can be deactivated
    % during the call. For future Matlab releases this way has to be
    % reevaluated.
    %
    % msgid:  'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame'
    %
    % msgstr: 'figure JavaFrame property will be obsoleted in a future
    %          release. For more information see <a
    %          href="http://www.mathworks.com/javaframe">the JavaFrame
    %          resource on the MathWorks web site</a>.'
    if exist('verLessThan', 'file') && ~verLessThan('matlab', '9.2')
        % at least Matlab R2017a
        warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        maximize(hFigureHandle);
        warning('on', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    else
        maximize(hFigureHandle);
    end
    
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
