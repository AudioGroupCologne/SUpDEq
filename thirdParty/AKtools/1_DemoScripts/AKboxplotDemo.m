% This script show exemplary usage of AKboxplot, a tool for plotting data
% by means of descriptive statistics such as mean oder median, confidence
% intervals, and quartile ranges. This script shows the basic usage, for
% more information see the help inside AKboxplot
%
% 1. default box plot
% 2. adjust labels, and axis
% 3. change plot layout
% 4. change the box layout and categories
% 5. highlight deviation from zero, and show additional information
% 6. output data
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
close all; clear; clc

% get ratings - these are fake ratings for
% 60 subjects
% 8  dependen variables (e.g. perceptual qualities)
% 3  test conditions (e.g. audio contents)
getRatings


%%  --------------------------------------------------- 1. default box plot

% use default values for plotting. This will show:
% - boxes specifying the inter quartile range
% - lines and whiskers showing the 95% percent confidence interval
% - dots showing outliers
AKboxplot(ratings)

%% --------------------------------------------- 2. adjust labels, and axis
% parameters are passed to AKboxplot as a struct

% these are the lables for the x-axis - we make something up here
p.label = {'Difference' 'Liking' 'Naturalness' 'Tone Color' 'Location' 'Externalization' 'Extension' 'Dynamic'};

% we also adjust the y-axis layout
p.yLim  = [-1.2 1.2];
p.yTick = -1:1;

AKboxplot(ratings, p)

%% -------------------------------------------------- 3. change plot layout
% there are lots of parameters to cater the plot look to your needs

% we will have a plot with two rows - e.g. to use in a double column paper
p.rows = {1:3; 4:8};
% this adjust the plot size in cm
p.figure = [15 15];

% NOTE: see the help of AKboxplot for instructions on how to adjust the
% positioning of the subplots - this is done using subtightplot

% we will adjust the size of the dots showing the outliers and the
% linewidth for the smaller plot size
p.outlierSize = 10;
p.boxLine     = [1 1 2];

AKboxplot(ratings, p)

%% -------------------------------- 4. change the box layout and categories
% if your data is not normally distributet you usally don't show the
% confidence intervalls. We will use the lines to show the range instead
p.box = {'range' false 'iqr'};

% we can also change the color of the boxes using rgb vectors or strings
p.boxColor  = {'k'; 'w'; 'b'};
p.lineColor = {'none'; 'b'; 'none'};

AKboxplot(ratings, p)

%% ------ 5. highlight deviation from zero, and show additional information
% Data can be higlighted to make a plot more readable. This is illustrated
% using parts of the ratings.

% this will show two boxes. The first shows the 10-90 percentile range and
% the second the inter quartile range.
p2.box = {'range', [90 10], 'iqr'};

% this will color boxes gray, if the inter quartile range overlaps with
% zero
p2.criticalBound = 0;
p2.criticalValue = 'iqr';
p2.boxColor = {'k', [.6 .6 .6]};

AKboxplot(ratings(:,:,1), p2);

%% --------------------------------------------------------- 6. output data

% AKboxplot returns
% - dataOut: the plotted values (e.g. confidence intervals) in a table or
%            struct
% - specs  : complete specifications used for plotting
% - handles: handles to the subplot axis
[dataOut, specs, handles] = AKboxplot(ratings, p);