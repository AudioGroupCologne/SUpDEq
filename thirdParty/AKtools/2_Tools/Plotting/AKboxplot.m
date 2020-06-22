% [dataOut, specs, handles] = AKboxplot(dataIn, specs)
% can be used for highly modifiably boxplots
%
% for example usage see AKboxplotDemo.m
%
% I N P U T
% dataIn - data to be plotted. Can be a vector, or two or three dimensional
%          matrix of size [L M N], where L is the number of observations
%          for each case (e.g. number of subjects in a listening test), M
%          is the number of observed variables (e.g. perceptual items), and
%          N is the number of test conditions (e.g. audio contents).
%          AKboxplot will show M*N bars, and will group the N bars for each
%          observed variable. NaN-values will be omitted.
% specs  - a struct that defines the plot style with the fields listed
%          below. All fields are optionally, and the default parameters are
%          used if a field does not exist inside specs
%
%  label        : a cellstring array of size {M 1} that contains the names
%                 of the ovserved variables. This is used for labeling the
%                 x-axis (default = empty cell array).
%  figure       : a two element vector containing the width and height of
%                 the figure in cm. If a figure exists, this is ignored and
%                 the current figure is used for plotting. By default a
%                 large figure is created using AKf().
%  rows         : a cell array of size {R 1} that determines how many rows
%                 (subplots) are used for plotting the data.
%                 e.g. {1:4; [5 6 7 9]} plots variables 1-4 in the first
%                 row, and variables 5, 6, 7, and 9 in the second. All
%                 numbers in rows must <= M (see dataIn).
%                 specs.rows can contain different number of elements in
%                 each row. By default no subplots are used an all data is
%                 plotted in a single row.
%  columns      : a two element vector sepcifying the number of columns in
%                 the subplot layout, and the current column in which the
%                 data is plotted. By default only one column is used.
%  center       : true, or false to specify if boxes are centered within
%                 each rows (defualt = true).
%  categories   : Additional labels that are shown inside the plot. This is
%                 a cell array of cell arrays where each row contains as
%                 many entries as the same row in specs.rows.
%                 e.g. {{'OVERALL' '' ''}; {'COLORATION' '' '' ''}}
%                 hold information for a boxplot with two rows, with 3
%                 variables plotted in the first, and 4 variables plotted
%                 in the second row (default is an emtpy cell array).
%  categoriesY  : y-position of the labels from specs.categories. Can be
%                 'top' (default), 'bottom', or a scalar that specifies the
%                 y-value. Text will be shown in italic font.
%  groupSpacing : space between groups of bars in case N>1 (see dataIn).
%                 Scalar between [1 inf[ where 1 results in the same
%                 spacing as within groups of bars, and 2 will double the
%                 spacing (default = 2).
%  box          : 1 by 3 or N by 3 cell array that specifies the layout of
%                 the bars. The value inside the first column will be
%                 plotted as a line, the value inside the second column as
%                 a narrow box, and the value inside the third column as a
%                 wide box. Possible values are:
%                 > 'range' : spans between the minimum and maximum
%                 > 'ci95'  : spans the 95% confidence interval. Note that
%                             'ci' can be followed by any number between 0
%                             and 100
%                 > 'std1'  : spans between 1 times the standard deviation.
%                             'std' can be followed by a scalar, e.g.
%                             'std2' spans 2 times the standard deviation.
%                 > 'custom': allows passing custom data. This is only
%                             allowed for the first column/entry. See input
%                             parameter 'range' for further specification.
%                 > 'iqr'   : spans the inter quartile range
%                 > [a b]   : spans between a% and b% percentile
%                 > false   : doese not show the element
%                 defualt is {'ci95' false 'iqr'}
%  boxTH        : scalar or 1 by N vector of integers. If a column in
%                 dataIn contains less than specs.boxTH valid values the
%                 raw data is plotted as dots instead of showing boxes.
%                 (default = false)
%  boxCenter    : specifies the box center
%                 > 'mean' (default), 'median', or false
%  boxWhiskers  : true (default) or false to show or hide whiskers
%  boxWidth     : 1 by 3 or N by 3 array with elements between [0 1] that
%                 specifies the width of the whiskers (first column) wide
%                 box (second column) and specs.boxCenter (third column).
%                 Defult is [.5 .5 1].
%  boxLine      : 1 by 3 or N by 3 array that defines the line width of the
%                 whiskers and range (first column), boxes (second column),
%                 and specs.boxCenter (third column). Default is [1.5 0 3].
%  boxOutlier   : true (default) or false to show or hide outliers.
%                 Outliers are points in dataIn that are outside the range
%                 given by the first column of specs.box
%  outlierSize  : marker size for plotting specs.boxOutlier (Default = 20).
%  range        : Upper and lower range given in [2 M N] matrix (see input
%                 parameter 'dataIn'). 
%  criticalBound: boxes and lines can have different colors, if the do not
%                 overlap with the range given by specs.criticalBound. This
%                 is either a scalar, a two element vector, or false
%                 (default).
%  criticalValue: string that specifies what values are checked against
%                 specs.criticalBound. specs.criticalValue can have the
%                 same vales as specs.box.
%  boxColor     : Specification of box color. This can either be a single
%                 element, or a [1 by 2] cell array, an a [N by 1] cell
%                 array, or a [N by 2] cell array. Each row holds the color
%                 for a condition of dataIn. The first column holds the
%                 colors that are used if specs.criticalValue is outside
%                 specs.criticalBound, and the second column the colors in
%                 case specs.criticalValue overlaps specs.criticalBound.
%                 Each element can be either be a string (e.g. 'k' for
%                 black) or a three element RGB vector. By default,
%                 different colors are picked for all conditions
%  lineColor    : specifies the line color. Format identical to
%                 specs.boxColor, except that 'none' can be passed for not
%                 plotting lines (default = 'none').
%  yLabel       : string specifying the y label (default = false)
%  yLim         : two element vector specifying the y-axis limits in
%                 ascending order. By default the minimum and maximum
%                 vlaues of dataIn are used
%  yTick        : vector containing the y-ticks in ascending order. By
%                 default specs.yLim is used and zero is added if within
%                 the range of specs.yLim.
%  yTickLabel   : vector or cell array of same size as specs.yTicks. By
%                 default specs.yTick dataIn is plotted into the first
%                 column of the subplot layout and an empty cell is used
%                 otherwise.
%  yLines       : y-coordinate of intermediate lines to be plottet (scalar,
%                 row vector, or column vector. By default a line is drawn
%                 at 0.
%  yDash        : false or true to specifies if grid lines are dashed. Pass
%                 a scalar or vector to specify the lengths of dashes and
%                 white spaces (ws) in between
%                 > [a]       : dashes and ws of a mm
%                 > [a b]     : dashes of a, and ws of b mm
%                 > [a b c d] : dash if a, ws of b, dash of c, and ws of d
%                               mm
%  gap          : subtightplot is used for the subplot layout. See doc
%                 subtightplot.m (default = [.05 0])
%  margH        : subtightplot is used for the subplot layout. See doc
%                 subtightplot.m (default = [.075 .05])
%  margW        : subtightplot is used for the subplot layout. See doc
%                 subtightplot.m (default = [.15 .05], if we plot inside
%                 the first column of a subplot, [.075 .05] elsewise)
%
% O U T P U T
% dataOut - specifies the ranges used for plotting the boxes for all
%           variables and conditions of dataIn
% specs   - specs complemented by default parameters
% handles - handles to the subplot axes
%
%
% 12/2016 - fabian.brinkmann@tu-berlin.de

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
function [dataOut, specs, handles] = AKboxplot(dataIn, specs)
%% ------------------------------------------------------- 1. data handling
if size(dataIn, 1)==1
    dataIn = dataIn';
end

if nargin==1 || ~isstruct(specs)
    specs = {};
end

% ------------------------------ set default parameter and check input data
if ~isfield(specs, 'label')
    specs.label = cell(size(dataIn, 2), 1);
end

if ~isfield(specs, 'rows')
    specs.rows = {1:size(dataIn, 2)};
end

if ~isfield(specs, 'columns')
    specs.columns = [1 1];
end


if ~isfield(specs, 'center')
    specs.center = true;
end

if ~isfield(specs, 'categories')
    for nn = 1:numel(specs.rows)
        specs.categories{nn,1} = cell(numel(specs.rows{nn}),1);
    end
end

if ~isfield(specs, 'categoriesY')
    specs.categoriesY = 'bottom';
end

if ~isfield(specs, 'groupSpacing')
    specs.groupSpacing = 2;
end

if ~isfield(specs, 'boxTH')
    specs.boxTH = zeros(1, size(dataIn,3));
elseif numel(specs.boxTH) < size(dataIn,3)
    specs.boxTH = repmat(specs.boxTH(1), [1 size(dataIn,3)]);
end

if ~isfield(specs, 'box')
    specs.box = {'ci95' false 'iqr'};
    specs.box = repmat(specs.box, [size(dataIn,3) 1]);
elseif size(specs.box,1) < size(dataIn,3)
    specs.box = repmat(specs.box(1,:), [size(dataIn,3) 1]);
end

if ~isfield(specs, 'boxCenter')
    specs.boxCenter = {'mean'};
    specs.boxCenter = repmat(specs.boxCenter, [1 size(dataIn,3)]);
elseif ischar(specs.boxCenter) || islogical(specs.boxCenter)
    specs.boxCenter = repmat({specs.boxCenter}, [1 size(dataIn,3)]);
end

if ~isfield(specs, 'boxWhiskers')
    specs.boxWhiskers = true(1, size(dataIn,3));
elseif numel(specs.boxWhiskers) < size(dataIn,3)
    specs.boxWhiskers = repmat(specs.boxWhiskers, [1 size(dataIn,3)]);
end

if ~isfield(specs, 'boxOutlier')
    specs.boxOutlier = true(1, size(dataIn,3));
elseif numel(specs.boxOutlier) < size(dataIn,3)
    specs.boxOutlier = repmat(specs.boxOutlier, [1 size(dataIn,3)]);
end

if ~isfield(specs, 'boxWidth')
    specs.boxWidth = repmat([.5 .5 1], [size(dataIn, 3) 1]);
elseif size(specs.boxWidth, 1) < size(dataIn,3)
    specs.boxWidth = repmat(specs.boxWidth(1,:), [size(dataIn,3) 1]);
end

if ~isfield(specs, 'boxLine')
    specs.boxLine = repmat([1.5 0 3], [size(dataIn, 3) 1]);
elseif size(specs.boxLine, 1) < size(dataIn,3)
    specs.boxLine = repmat(specs.boxLine(1,:), [size(dataIn,3) 1]);
end

if ~isfield(specs, 'range')
    specs.range = false;
end

if ~isfield(specs, 'outlierSize')
    specs.outlierSize = 20*ones(1, size(dataIn,3));
elseif numel(specs.outlierSize) < size(dataIn,3)
    specs.outlierSize = repmat(specs.outlierSize(1), [1 size(dataIn,3)]);
end

if ~isfield(specs, 'criticalBound')
    specs.criticalBound = false;
end

if ~isfield(specs, 'criticalValue')
    specs.criticalValue = false;
else
    if ~islogical(specs.criticalValue) && islogical(specs.criticalBound)
        specs.criticalBound = 0;
    end
end

if ~isfield(specs, 'boxColor')
    % possible colors
    c = 'kbrypgldmc';
    % create cell array with colrs
    specs.boxColor = cell(size(dataIn,3), 1);
    for nn = 1:size(dataIn,3)
        pick = mod(nn,10);
        if ~pick
            pick = 10;
        end
        specs.boxColor{nn} = AKcolors(c(pick));
    end
    
    specs.boxColor = repmat(specs.boxColor, [1 2]);
    
    clear c nn pick
else
    if ~iscell(specs.boxColor)
        if ischar(specs.boxColor)
            specs.boxColor = AKcolors(specs.boxColor);
        end
        
        specs.boxColor = repmat({specs.boxColor}, [size(dataIn,3) 2]);
        
    else
        if size(specs.boxColor, 2) == 1
            specs.boxColor = repmat(specs.boxColor, [1 2]);
        end
        if size(specs.boxColor, 1) < size(dataIn,3)
            specs.boxColor = repmat(specs.boxColor(1,:), [size(dataIn,3) 1]);
        end
    end
end

if ~isfield(specs, 'lineColor')
    specs.lineColor = specs.boxColor;
else
    if ~iscell(specs.lineColor)
        if ischar(specs.lineColor)
            if ~strcmpi(specs.lineColor, 'none')
                specs.lineColor = AKcolors(specs.lineColor);
            end
        end
        
        specs.lineColor = repmat({specs.lineColor}, [size(dataIn,3) 2]);
        
    else
        if size(specs.lineColor, 2) == 1
            specs.lineColor = repmat(specs.lineColor, [1 2]);
        end
        if size(specs.lineColor, 1) < size(dataIn,3)
            specs.lineColor = repmat(specs.lineColor(1,:), [size(dataIn,3) 1]);
        end
    end
end

if ~isfield(specs, 'yLabel')
    specs.yLabel = false;
end

if ~isfield(specs, 'yLim')
    specs.yLim = [floor(min(dataIn(:))) ceil(max(dataIn(:)))];
end

if ~isfield(specs, 'yTick')
    if specs.yLim(1)<0 && specs.yLim(2) > 0
        specs.yTick = [specs.yLim(1) 0 specs.yLim(2)];
    else
        specs.yTick  = specs.yLim;
    end
else
    specs.yTick = sort(specs.yTick);
end

if ~isfield(specs, 'yTickLabel')
    if specs.columns(2) == 1
        specs.yTickLabel = specs.yTick;
    else
        specs.yTickLabel = '';
    end
end

if ~isfield(specs, 'yLines')
    specs.yLines = specs.yTick;
end

if ~isfield(specs, 'yDash')
    specs.yDash = [2 2 2 2];
else
    if islogical(specs.yDash) && specs.yDash
        specs.yDash = [2 2 2 2];
    elseif numel(specs.yDash) ~= 4
        specs.yDash = repmat(specs.yDash, [1 4/numel(specs.yDash)]);
    end
end

if ~isfield(specs, 'gap')
    specs.gap = [.05 0];
end

if ~isfield(specs, 'margH')
    specs.margH = [.075 .05];
end

if ~isfield(specs, 'margW')
    if specs.yLabel
        specs.margW = [.15 .05];
    else
        specs.margW = [.075 .05];
    end
end

%% --------------------------------------------------- 2. helping variables

% maximum number of items in one row
sub_n = zeros(numel(specs.rows),1);
for nn = 1:numel(specs.rows)
    sub_n(nn) = numel(specs.rows{nn});
end
subN = max(sub_n);

% convert colors from strings to RGB
for nn = 1:size(specs.boxColor,1)
    for mm = 1:size(specs.boxColor,2)
        if ischar(specs.boxColor{nn,mm})
            specs.boxColor{nn,mm} = AKcolors(specs.boxColor{nn,mm});
        end
        if ischar(specs.lineColor{nn,mm})
            if ~strcmpi(specs.lineColor{nn,mm}, 'none')
                specs.lineColor{nn,mm} = AKcolors(specs.lineColor{nn,mm});
            end
        end
        
        if (all(specs.boxColor{nn,mm} == 1) && all(specs.lineColor{nn,mm} == 1)) || ...
           (all(specs.boxColor{nn,mm} == 1) && strcmpi(specs.lineColor{nn,mm}, 'none'))
            specs.lineColor{nn,mm} = [0 0 0];
        end
    end
end

clear nn

%% ------------------------------------------------------------- 3. boxplot

% allocate output data
writeOut  = 0;
desc1     = nan(sum(sub_n)*size(dataIn,3), 2);
desc2     = desc1;
desc3     = desc1;
center    = nan(sum(sub_n)*size(dataIn,3), 1);
itemNames = cell(sum(sub_n)*size(dataIn,3), 1);
condition = nan(sum(sub_n)*size(dataIn,3), 1);
handles   = zeros(numel(specs.rows), 1);

% make figure if none exist
if isempty(findall(0,'Type','Figure'))
    if isfield(specs, 'figure')
        AKf(specs.figure(1),specs.figure(2))
    else
        AKf
    end
else
    specs.figure = 'parent';
end

% loop over subplots
for nn = 1:numel(specs.rows)
    
    % ------------------------------------- get item IDs in current subplot
    items = specs.rows{nn};
    
    % ------------------------------------------------ set subplot and axis
    subtightplot(numel(specs.rows), specs.columns(1), specs.columns(1)*(nn-1)+specs.columns(2), specs.gap, specs.margH, specs.margW)
    handles(nn) = gca;
    
    % maximum x-value
    xMax = subN*size(dataIn,3)+(subN-1)*specs.groupSpacing+.5;
    axis([0.5 xMax min(specs.yLim) max(specs.yLim)])
    hold on
    box on
    
    % set xAxis ticks
    groupStart   = size(dataIn,3)/2+0.5;                                                % x-value for first xTick
    groupSpacing = size(dataIn,3)+specs.groupSpacing;                                   % spacing between xTicks
    groupTicks   = groupStart:groupSpacing:groupStart+(numel(items)-1)*groupSpacing;    % acutal xTicks
    % change xTicks if we have less than nSub items in the current subplot and want to center
    if specs.center && numel(items)<subN
        freeSpace = (subN-1)*groupSpacing - (groupStart+(numel(items)-1)*groupSpacing);
        freeSpace = freeSpace/2;
    else
        freeSpace = 0;
    end
    groupTicks = groupTicks+freeSpace;
    
    % set the values
    set(gca, 'xTick', groupTicks, ...
             'xTickLabel', specs.label(items), ...
             'yTick', specs.yTick, ...
             'yTickLabel', specs.yTickLabel, ...
             'TickLength', [0 0], ...
             'fontWeight', 'bold')
    
    % ylabel
    if any(specs.yLabel)
        ylabel(specs.yLabel)
    end
    
    % draw tick lines
    grid off
    for mm = 1:numel(specs.yLines)
        if all(specs.yLines(mm)~=specs.yLim)
            if specs.yDash
                dashline([0 xMax], [specs.yLines(mm) specs.yLines(mm)], specs.yDash(1), specs.yDash(2), specs.yDash(3), specs.yDash(4), 'k')
            else
                plot([0 xMax], [specs.yLines(mm) specs.yLines(mm)], 'k')
            end
        end
    end
    
    % --------------------------------------------------------- actual plot
    % loop over items
    for ii = 1:numel(items)        
        
        % ----------------------------------------------- categorial labels
        if cell2mat(specs.categories{nn}(ii))
            if strcmpi(specs.categoriesY, 'bottom')
                text(groupTicks(ii), specs.yLim(1), cell2mat(specs.categories{nn}(ii)), 'fontAngle', 'italic', 'horizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
            elseif strcmpi(specs.categoriesY, 'top')
                text(groupTicks(ii), specs.yLim(2), cell2mat(specs.categories{nn}(ii)), 'fontAngle', 'italic', 'horizontalAlignment', 'center', 'VerticalAlignment', 'top')
            else
                text(groupTicks(ii), specs.categoriesY, cell2mat(specs.categories{nn}(ii)), 'fontAngle', 'italic', 'horizontalAlignment', 'center')
            end
        end
        
        % loop over conditions
        for cc = 1:size(dataIn,3)
            
            writeOut = writeOut + 1;
            itemNames{writeOut} = specs.label{items(ii)};
            condition(writeOut) = cc;
            
            % current ratings
            ratings = dataIn(:, items(ii), cc);
            % remove unwanted values
            ratings = ratings(~isnan(ratings) & ~isinf(ratings));
            
            % current x-value for plotting
            x = (ii-1)*size(dataIn,3)+cc + (ii-1)*specs.groupSpacing + freeSpace;
            
            if numel(ratings) >= specs.boxTH(cc) % <- plot boxes
                
                % ----- get the value for central tendency -----
                if any(specs.boxCenter{cc})
                    if strcmpi(specs.boxCenter{cc}, 'median')
                        yCenter  = prctile(ratings, 50);
                    elseif strcmpi(specs.boxCenter{cc}, 'mean')
                        yCenter  = mean(ratings);
                    end
                    % save to output
                    center(writeOut) = yCenter;
                else
                    yCenter  = prctile(ratings, 50);
                end
                
                % ----- get the color -----
                if any(specs.criticalValue)
                    
                    % get the range for comparing it to the critical value
                    if isnumeric(specs.criticalValue)
                        cr = prctile(ratings, specs.criticalValue);
                    elseif strcmpi(specs.criticalValue, 'iqr')
                        cr = prctile(ratings, [75 25]);
                    elseif strcmpi(specs.criticalValue, 'range')
                        cr = prctile(ratings, [100 0]);
                    elseif strcmpi(specs.criticalValue(1:2), 'ci')
                        if isempty( specs.criticalValue(3:end) )
                            cr = AKci(ratings, 95);
                        else
                            cr = AKci(ratings, str2double(specs.criticalValue(3:end)));
                        end
                    elseif ~isempty( strfind( lower(specs.criticalValue), 'std') )
                        if isempty( specs.criticalValue(4:end) )
                            cr = yCenter + [-std(ratings) std(ratings)];
                        else
                            cr = yCenter + [-std(ratings) std(ratings)]*str2double( specs.criticalValue(4:end) );
                        end
                    elseif strcmpi(specs.criticalValue, 'custom')
                        cr = specs.range(:,ii,cc)';
                    else
                        error('AKboxplot:Input', [specs.criticalValue ' is not a valid value for ''criticalValue'''])
                    end
                    
                    % check if we overlap the critical value or critical range
                    if numel(specs.criticalBound) == 1
                        if all(cr>max(specs.criticalBound)) || all(cr<min(specs.criticalBound))
                            colorPick = 1;
                        else
                            colorPick = 2;
                        end
                    else
                        if ~( any(cr<max(specs.criticalBound)) && any(cr>min(specs.criticalBound)) )
                            colorPick = 1;
                        else
                            colorPick = 2;
                        end
                    end
                    
                else
                    colorPick = 1;
                end
                
                % pick the color and alpha
                c  = specs.boxColor{cc,colorPick};
                lc = specs.lineColor{cc,colorPick};
                
                % use different color for plotting the range and whiskers
                % if the line color is none or white
                if strcmpi(lc, 'none') || all(lc == 1)
                    if all(c==1)
                        lcRange = [0 0 0];
                    else
                        lcRange = c;
                    end
                else
                    lcRange = lc;
                end
                
                % ----- get the y-values for plotting lines and boxes -----
                for mm = 1:3
                    % get descriptives
                    if any(specs.box{cc,mm})
                        if isnumeric(specs.box{cc,mm})
                            y = prctile(ratings, specs.box{cc,mm});
                        elseif strcmpi(specs.box{cc,mm}, 'IQR')
                            y = prctile(ratings, [75 25]);
                        elseif strcmpi(specs.box{cc,mm}, 'range')
                            y = prctile(ratings, [100 0]);
                        elseif ~isempty( strfind( lower(specs.box{cc,mm}), 'std') )
                            if isempty( specs.box{cc,mm}(4:end) )
                                y = yCenter + [-std(ratings) std(ratings)];
                            else
                                y = yCenter + [-std(ratings) std(ratings)]*str2double( specs.box{cc,mm}(4:end) );
                            end
                        elseif strcmpi(specs.box{cc,mm}, 'custom')
                            y = specs.range(:,ii,cc)';
                        elseif strcmpi(specs.box{cc,mm}(1:2), 'ci')
                            if isempty( specs.box{cc,mm}(3:end) )
                                y = AKci(ratings, 95);
                            else
                                y = AKci(ratings, str2double(specs.box{cc,mm}(3:end)));
                            end
                        else
                            error('AKboxplot:Input', [specs.box{cc,mm} ' is not a valid value for ''box'''])
                        end
                    else
                        y = false;
                    end
                    
                    % assign to correct variable
                    if mm == 1
                        yLine = y;
                        desc1(writeOut,:) = y;
                    elseif mm == 2
                        yBox1 = y; %#ok<NASGU>
                        desc2(writeOut,:) = y;
                    else
                        yBox2 = y; %#ok<NASGU>
                        desc3(writeOut,:) = y;
                    end
                end
                
                % ----- plot the data -----
                % plot central tendency
                if any(specs.boxCenter{cc})
                    plot([-specs.boxWidth(cc,3) specs.boxWidth(cc,3)]/2+x, [yCenter yCenter], 'color', lcRange(1:3), 'lineWidth', specs.boxLine(cc,3))
                end
                % plot line
                if any(specs.box{cc,1})
                    plot([x x], yLine, 'color', lcRange(1:3), 'lineWidth', specs.boxLine(cc,1))
                end
                if specs.boxWhiskers(cc)
                    plot([-specs.boxWidth(cc,1) specs.boxWidth(cc,1)]./4 + x, [yLine(1) yLine(1)], 'color', lcRange(1:3), 'lineWidth', specs.boxLine(cc,1))
                    plot([-specs.boxWidth(cc,1) specs.boxWidth(cc,1)]./4 + x, [yLine(2) yLine(2)], 'color', lcRange(1:3), 'lineWidth', specs.boxLine(cc,1))
                end
                % plot boxes
                for bb = 1:2
                    if any(specs.box{cc,bb+1})
                        % get the y-coordinates and box width
                        yData = eval(['yBox' num2str(bb)]);
                        if bb == 1
                            scale = 4;
                        else
                            scale = 2;
                        end
                        % plot the boxes with, or without lines
                        if specs.boxLine(cc,2)
                            patch([-specs.boxWidth(cc,2) -specs.boxWidth(cc,2) specs.boxWidth(cc,2) specs.boxWidth(cc,2)]/scale+x, [min(yData) max(yData) max(yData) min(yData)], c(1:3), ...
                                  'EdgeColor', lc(1:3), 'lineWidth', specs.boxLine(cc,2))
                        else
                            patch([-specs.boxWidth(cc,2) -specs.boxWidth(cc,2) specs.boxWidth(cc,2) specs.boxWidth(cc,2)]/scale+x, [min(yData) max(yData) max(yData) min(yData)], c(1:3), ...
                                  'EdgeColor', lc(1:3))
                        end
                    end
                end
                % plot outlier
                if specs.boxOutlier(cc)
                    yOutlier = ratings(ratings<min(yLine) | ratings>max(yLine));
                    if ~isempty(yOutlier)
                        plot(repmat(x, [1, numel(yOutlier)]), yOutlier, '.', 'color', lcRange(1:3), 'markerSize', specs.outlierSize(cc))
                    end
                end
                
            else % <- plot data as points only
                
                c = specs.boxColor{cc,1};
                plot(repmat(x, [1 numel(ratings)]), ratings, '.', 'color', c, 'markerSize', specs.outlierSize(cc))
                
            end
            
        end
    end
    
end

%% -------------------------------------------------- 4. format output data

if nargout
    % ---------------------------------------------- names of column in dataOut
    names = cell(4,1);
    for nn = 1:3
        if any(specs.box{1,nn})
            if isnumeric(specs.box{1,nn})
                names{nn} = ['prct_' num2str(specs.box{1,nn}(1)) '_' num2str(specs.box{1,nn}(2))];
            else
                names{nn} = specs.box{1,nn};
            end
        else
            names{nn} = ['none_' num2str(nn)];
        end
    end
    
    if any(specs.boxCenter{1})
        names{4} = specs.boxCenter{1};
    else
        names{4} = 'none';
    end
    
    names = {'Item' 'Condition' names{4} names{1} names{2} names{3}};
    
    % ------------------------------ write table or cell array with the results
    if exist('table.m', 'file')
        % create table
        dataOut = table(itemNames, condition, center, desc1, desc2, desc3, ...
            'variableNames', names);
        % remove empty columns
        if ~any(specs.boxCenter{1})
            dataOut.none = [];
        end
        if ~any(specs.box{1,1})
            dataOut.none_1 = [];
        end
        if ~any(specs.box{1,2})
            dataOut.none_2 = [];
        end
        if ~any(specs.box{1,3})
            dataOut.none_3 = [];
        end
    else
        dataOut.itemNames = itemNames;
        dataOut.condition = condition;
        dataOut.(names{3})    = center;
        dataOut.(names{4})     = desc1;
        dataOut.(names{5})     = desc2;
        dataOut.(names{6})     = desc3;
    end
end