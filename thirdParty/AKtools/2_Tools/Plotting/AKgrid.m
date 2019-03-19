function AKgrid(h, gridColor, gridLinewidth, dashGrid, drawOnTop)
% AKgrid(h, gridColor, gridLinewidth, dashGrid, drawOnTop)
%
% completely redraws grid of a figure, because the Matlab grids sometimes
% look ugly in plots. All input arguments are optional
%
% See AKplotDemo.m for examples
%
% I N P U T (default value):
% h (gca)                - handle to the axis the grid will be drawn in
% gridColor ([.8 .8 .8]) - grid line color. Three element RGB vector, or
%                          string (e.g. 'k')
% gridLinewidth (.5)     - grid line width
% dashGrid (false)       - draw dashed grid. False or four element vector
%                          [dashLength gap dashLength gap] with values in
%                          millimeter
% drawOnTop (false)      - draws the grid below the plotted data (default)
%                          or on top of it
%
% 09/2013 - fabian.brinkmann@tu-berlin.de
% 10/2016 - fabian.brinkmann@tu-berlin.de (added dashGrid, removed custom box)

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

% set default parameters
if ~exist('h', 'var')
    h = gca;
end
if ~exist('gridColor', 'var')
    gridColor = [.8 .8 .8];
end
if ~exist('gridLinewidth', 'var')
    gridLinewidth = .5;
end
if ~exist('dashGrid', 'var')
    dashGrid = false;
end
if ~exist('drawOnTop', 'var')
    drawOnTop = false;
end

% flip line order (it will be flipped later again...)
if ~drawOnTop
    set(h,'children',flipud(get(gca,'children')))
end

% get axes properties
%   for some weird reason matlab sometimes gets the wrong ticks if it does
%   not wait for his soul to rest - may peace be with you my friend
pause(.1)   
xTick = get(h, 'xTick');
xLim  = get(h, 'xLim');
yTick = get(h, 'yTick');
yLim  = get(h, 'yLim');

% if a log axis has 0 as smalles value, the lines will not be plotted
scale = get(h, 'xScale');
if strcmpi(scale, 'log')
    xTick(xTick == 0) = .1;
    xLim(xLim == 0) = .1;
end

scale = get(h, 'yScale');
if strcmpi(scale, 'log')
    yTick(yTick == 0) = .1;
    yLim(yLim == 0) = .1;
end

% remove grid and box
grid off

% redraw grid
hold on
for n = 1:length(xTick)
    if xTick(n) > xLim(1) && xTick(n) < xLim(2)
        if any(dashGrid)
            dashline(h, [xTick(n) xTick(n)], [yLim(1) yLim(2)], dashGrid(1), dashGrid(2), dashGrid(3), dashGrid(4), 'color', gridColor, 'Linewidth', gridLinewidth)
        else
            plot(h, [xTick(n) xTick(n)], [yLim(1) yLim(2)], 'color', gridColor, 'Linewidth', gridLinewidth)
        end
    end
end

for n = 1:length(yTick)
    if yTick(n) > yLim(1) && yTick(n) < yLim(2)
        if any(dashGrid)
            dashline(h, [xLim(1) xLim(2)], [yTick(n) yTick(n)], dashGrid(1), dashGrid(2), dashGrid(3), dashGrid(4), 'color', gridColor, 'Linewidth', gridLinewidth)
        else
            plot(h, [xLim(1) xLim(2)], [yTick(n) yTick(n)], 'color', gridColor, 'Linewidth', gridLinewidth)
        end
    end
end

% flip line order to move the grid to the background
if ~drawOnTop
    set(h,'children',flipud(get(gca,'children')))
end

% we don't need tick marks if we have a grid
set(h, 'Ticklength', [0 0])

% make a box in case there is non
box on

end






%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% downloaded from http://de.mathworks.com/matlabcentral/fileexchange/1892-dashline
% no licence provided :(
function [h1, h2]=dashline(h, xdata, ydata, dash1, gap1, dash2, gap2, varargin)
%DASHLINE  Function to produce accurate dotted and dashed lines
%   DASHLINE(XDATA, YDATA, DASH1, GAP1, DASH2, GAP2, ...) plots 
%   the data XDATA and YDATA usig a two dash linestyle. The 
%   function allows the lengths of the dashes to be specified. 
%   DASH1 is the length of the first dash; GAP1 is the length of 
%   the first gap; DASH2 is the length of the second  dash and 
%   GAP2 is the length of the second gap, with all lengths being 
%   in millimeters. The two dash pattern is repeated along the 
%   length of the line. The routine is designed to allow greater
%   control over the linestyle than is available using PLOT.
%
%   Additional inputs may be given to specify line properties. 
%   These are passed through to the plot function. 
%
%   DASH1 or DASH2 may be strings specifying plot symbols. A 
%   dashlength of 0 will cause a dot to be plotted. If GAP1 and 
%   GAP2 are 0 then the arguments are simply passed through to a 
%   plot command
%
%   For example: 
%   clf;
%   dashline([1:10],rand(1,10),1,1,1,1) 
%      %produces a dashed line with 1mm dashes and 1mm gaps
%   hold on;
%   dashline([1:10],rand(1,10),4,2,'o',2,'markerfacecolor','y',...
%      'color','k') 
%      %produces a line with black dashes and yellow centred circles
%
%   [H1, H2]=DASHLINE(...) Outputs are handles to the plotted 
%      lines and markers
%
%   DASHLINE works by calculating the distance along the line to be 
%   plotted, using the existing axes position. If HOLD is OFF when 
%   DASHLINE is called, then the axis limits are set automatically; 
%   otherwise the existing limits are used. DASHLINE then 
%   interpolates XDATA and YDATA at the positions of the start and 
%   end of the dashes, using NaN's to lift the pen. Rescaling the 
%   axes after DASHLINE has been called, or changing the limits, or 
%   changing from log to linear plots, will mean that the dashes no 
%   longer have the specified lengths. This routine does not calculate 
%   distances correctly if the axes DataAspectRatio or 
%   PlotBoxAspectRatio have been manually specified. The lines will 
%   not display properly in a call to LEGEND.
%
%   Kluged together by Edward Abraham (e.abraham@niwa.cri.nz) 
%   27-Sept-2002 Minor change made so that it works with Matlab 6.5
%   27-June-2002 Original version




%Check shape of input data, and put into column form ..
[si, sj]=size(xdata);
[siy, sjy]=size(ydata);
if (si*sj ~= siy*sjy || ((si~=1 && sj~=1) && (si~=siy || sj~=sjy)))
    error('Vectors must be the same lengths.')
end
if ( si==1 || siy==1 )
    xdata=xdata(:);
    ydata=ydata(:);
    si=sj;
    sj=1;
end

%Check inputs for the length of the dashes and gaps are sensible ...
if ischar(gap1) || ischar(gap2) || ~isreal(gap1) || ~isreal(gap2) || ~isfinite(gap1) || ~isfinite(gap2) || gap1<0 || gap2<0
    error('Gaps must be positive lengths')    
end
if (~ischar(dash1) && ( ~isreal(dash1) || ~isfinite(dash1) || dash1<0)) || (~ischar(dash2) && ( ~isreal(dash2) || ~isfinite(dash2) || dash1<0))
    error('Dashes must either be positive lengths or plot strings.')
end   
% If there are no gaps between the dashes then pass data straight through to a plot comand
if (gap1==0 && gap2==0)
    p=plot(xdata,ydata,'-',varargin{:});
    if (nargout>0)
        h1=p;
        h2=[];
    end
    return
end

%Check inputs for string specification of the dashes (dashs of zero length are plotted as dots) ...
if ischar(dash1)
    Marker1=dash1;
    dash1=0;
elseif dash1==0
    Marker1='.';
end
if ischar(dash2)
    Marker2=dash2;
    dash2=0;
elseif dash2==0
    Marker2='.';
end


%Get Axes properties ...
AxesUnits=get(h,'units');
IsXLog=strcmp(get(h,'XScale'),'log');
IsYLog=strcmp(get(h,'YScale'),'log');
IsHold=ishold;
set(h,'units','centimeters');
Position=get(h,'Position');

try
    % If hold is off then determine axes limits by initially plotting the data ..
    if ~IsHold
        cla
        p=plot(xdata, ydata);
        if IsXLog
            set(h,'Xscale','log');
        end
        if IsYLog
            set(h,'Yscale','log');
        end
    end
    hold on;
    XLim=get(h,'Xlim');
    YLim=get(h,'Ylim');
    if ~IsHold
        delete(p)
    end

    % Try to correct for the annoying fact that in Log mode the axes limits are not always the axes limits!
    if IsXLog && XLim(1)==0
        XTick=get(h,'xtick');
        XLim(1)=XTick(1);
    end
    if IsYLog && YLim(1)==0
        YTick=get(h,'ytick');
        YLim(1)=YTick(1);
    end


    %Work out position of datapoints...
    if ~IsXLog
        xpos=(xdata-XLim(1))/(XLim(2)-XLim(1))*Position(3);
    else
        xpos=(log10(xdata)-log10(XLim(1)))/(log10(XLim(2))-log10(XLim(1)))*Position(3);
    end

    if ~IsYLog
        ypos=(ydata-YLim(1))/(YLim(2)-YLim(1))*Position(4);
    else
        ypos=(log10(ydata)-log10(YLim(1)))/(log10(YLim(2))-log10(YLim(1)))*Position(4);
    end

    handles=NaN*ones(sj,3);
    %Process each column of data ...
    for i=1:sj
        xposi=xpos(:,i);
        yposi=ypos(:,i);
        xdatai=xdata(:,i);
        ydatai=ydata(:,i);
        f=find(~isreal(xposi) | isinf(xposi)  | isnan(xposi) | ~isreal(yposi) | isinf(yposi)  | isnan(yposi));
        if ~isempty(f)
            xposi(f)=[];
            yposi(f)=[];
            xdatai(f)=[];
            ydatai(f)=[];
        end
        %Calculate distance from the start of the line (in mm) ...
        dist=[0;cumsum(sqrt(diff(xposi).^2 + diff(yposi).^2))*10];

        start1=0:dash1+gap1+dash2+gap2:dist(end);
        dashes=zeros(6*length(start1),1);
        dashes(1:6:end)=start1;
        dashes(2:6:end)=start1+dash1;
        dashes(3:6:end)=NaN;
        dashes(4:6:end)=start1+dash1+gap1;
        dashes(5:6:end)=start1+dash1+gap1+dash2;
        dashes(6:6:end)=NaN;
    
        xdash=NaN*zeros(length(dashes)+length(xdata),1);
        ydash=NaN*zeros(length(dashes)+length(xdata),1);
        %Straight dashes ...
        if ~IsXLog
            xdash(1:length(dashes))=interp1(dist, xdatai, dashes);
        else
            xdash(1:length(dashes))=10.^interp1(dist, log10(xdatai), dashes);
        end 
        if ~IsYLog
            ydash(1:length(dashes))=interp1(dist, ydatai, dashes);
        else
            ydash(1:length(dashes))=10.^interp1(dist, log10(ydatai), dashes);
        end 
        % Get data for markers ...
        if (dash1==0)
            xdot1=xdash(1:6:end);
            ydot1=ydash(1:6:end);
            dashes(1:6:end)=NaN;
            dashes(2:6:end)=NaN;
            xdash(1:6:end)=NaN;
            xdash(2:6:end)=NaN;
            ydash(1:6:end)=NaN;
            ydash(2:6:end)=NaN;
        end
        if (dash2==0)
            xdot2=xdash(4:6:end);
            ydot2=ydash(4:6:end);
            dashes(4:6:end)=NaN;
            dashes(5:6:end)=NaN;
            xdash(4:6:end)=NaN;
            xdash(5:6:end)=NaN;
            ydash(4:6:end)=NaN;
            ydash(5:6:end)=NaN;
        end
    
        %Insert data points that fall within dashes (allows dashes to curve...)
        count=0;
        xlen=length(xdash);
        dashstart=zeros(2*length(start1),1);
        dashstart(1:2:end)=1:6:length(dashes);
        dashstart(2:2:end)=4:6:length(dashes);
        for j=1:length(dashstart)
            f=find(dist>dashes(dashstart(j)) & dist<dashes(dashstart(j)+1));
            if ~isempty(f)
                xdash(dashstart(j)+count+length(f)+1:xlen+length(f))=xdash(dashstart(j)+count+1:xlen);
                xdash(dashstart(j)+count+1:dashstart(j)+count+length(f))=xdatai(f);
                ydash(dashstart(j)+count+length(f)+1:xlen+length(f))=ydash(dashstart(j)+count+1:xlen);
                ydash(dashstart(j)+count+1:dashstart(j)+count+length(f))=ydatai(f);
                xlen=xlen+length(f);
                count=count+length(f);
            end
        end
    
    
        %Plot line and markers ...
        if (~(dash1==0 && dash2==0))
            handles(i,1)=plot(xdash, ydash, varargin{:});
        end
        if (dash1==0)
            handles(i,2)=plot(xdot1, ydot1,Marker1,varargin{:});
        end
        if (dash2==0)
            handles(i,3)=plot(xdot2, ydot2,Marker2,varargin{:});
        end
    end

    % Only return two handles (h2 is empty if there was only one call to plot for each column of data) ...
    if (nargout>0)
        handles(isnan(handles))=[];
        h1=handles(:,1);
        [~, hj]=size(handles);
        if (hj>1)
            h2=handles(:,2);
        else
            h2=[];
        end
    end

    %Restore Axes properties ...
    set(h,'Units',AxesUnits,'XLim',XLim,'YLim',YLim);
    if ~IsHold
        hold off
    end
catch ME
    %Restore Axes properties ...
    set(h,'Units',AxesUnits,'XLim',XLim,'YLim',YLim);
    if ~IsHold
        hold off
    end
    error(ME.stack)
end
end
