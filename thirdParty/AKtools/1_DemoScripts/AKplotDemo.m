% Demo script for the AK plotting tool AKp.m
%
% 1. 2D plots
% 2. 3D plots
% 3. spherical plots
% 4. polar plots
% 5. plots of multiple types
% 6. headphone impulse response (HpIR) plot
%
% v1.0 2012/05 fabian.brinkmann@tu-berlin.de

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

% get data for plotting (we use HRIRs for this)

% spatial sampling grid defined by azimuth and elevation (change false to
% true if you want to see the grid)
g = AKgreatCircleGrid(90:-10:-90, 10, 90, false);

% get left, and right ear HRIRs from spherical harmonics coefficients
AKdependencies('FABIAN')
[l, r] = AKhrirInterpolation(g(:,1), g(:,2), 0, 'measured_sh');

%% ------------------------------------------------------------ 1. 2D plots
%  AKp.m can be used for lots of different plots. In most cases it expects
%  impulse responses as input.
%  You can pass a string argument to change the plot type
%  't2d'  - time curve
%  'et2d' - log. time signal
%  'ed2d' - log. energy decay curve
%  'm2d'  - magnitude spectrum
%  'p2d'  - phase spectrum
%  'pu2d' - unwrapped phase spectrum
%  'gd2d' - group delay
%  'sr2d' - step response
%
%  you can play around with the next line to check out the different plot
%  types. This plots the HRIR for a source to the left of the listener

AKp([l(:,227) r(:,227)], 't2d')

%% you can quickly adjust the plot and make the grid look nice
%  (which already is the case with newer Matlab versions...)

AKf
AKp([l(:,227) r(:,227)], 'm2d', 'N', 44100/10, 'dr', [-20 20], 'x', [10 20000])
AKgrid

%% frequency domain plots can be smoothed

AKf
AKp([l(:,227) r(:,227)], 'm2d', 'frac', 3)

%% you can also get nice dashed lines for publications
%  thanks to Edward Abraham's dashline.m

AKf
AKp([l(:,227) r(:,227)], 'm2d', 'frac', 3, 'dash', [4 2])

%% There are some more plot types, that we did not check out so far
%  'toa2d' - broad band time of arrival (estimated from onset detection)
%  'itd2d' - broad band interaural time differences in HRIRs
%            (estimated from TOA)
%  'ild2d' - borad band interural level differences

% to show this we take only the horizontal plabe HRIRs
% (insert true instead of false to see the selected HRIRs)
id = AKsubGrid(g, 'transverse', 0, .1, false);
h.l = l(:,id);
h.r = r(:,id);


AKf(25,15) % <- takes care of your figure margins for saving to pdf
subplot(1,3,1)
AKp(h.l, 'toa2d', 'x', g(id,1), 'c', 'b')
AKp(h.r, 'toa2d', 'x', g(id,1), 'c', 'r')
subplot(1,3,2)
AKp(h, 'itd2d', 'x', g(id,1))
subplot(1,3,3)
AKp(h, 'ild2d', 'x', g(id,1))

%% The figure borders can be removed for plotting
AKtightenFigure


%% ------------------------------------------------------------ 2. 3D plots
% These plot types display lots of impulse responses and use color for
% showing the information you want to see
% 't3d', 'et3d', 'ed3d', 'sr3d', 'm3d', 'p3d', 'pu3d', 'g3d'

% once again we take our horizontal plane HRIRs
id = AKsubGrid(g, 'transverse', 0);

AKf
AKp(l(:,id), 'm3d', 'y', g(id,1))

%% you can change orientation of the plot and the colormap amongst others
%  see AKcolormaps for a documentation of color maps

AKf
AKp(l(:,id), 'm3d', 'y', g(id,1), 'hp_view', 'top_h', 'cm', 'AKwr')

%%  call AKcolormaps() for a list of brewer colormaps included in AKtools
%   (http://colorbrewer.org/).

AKcolormaps

%% you can quickly change the resolution of your colormap which comes
%  in handy from time to time

AKf
AKp(l(:,id), 'm3d', 'y', g(id,1), 'dr', [-20 20], 'cr', 2.5)  


%% ----------------------------------------------------- 3. spherical plots
%  These plot types can be used to plot data that is available on spherical
%  sampling grids, just like our HRIRs. They are available for
%  'm'   - magnitude spectra
%  'p'   - phase spectra
%  'pu'  - unwrapped phase spectra
%  'itd' - interaural time differences
%  'ild' - interaural level differences
%  'toa' - time of arrival
%
%  there are different display types of spherical plots
%  1: Balloon with fixed radius and color information
%  2: Balloon with variable radius and color information
%  3: Balloon variable radius and color decoding the phase
%  4: Balloon with variable radius but fixed color
%  5: Planar plot

% lets see how they look
AKf
AKp(l, 'm1', 'g', g, 'sph_f', 8000)

% The colored axis show the orientation of the data set:
% x-axis is red, point marker denotes poitive x
%   positive x points to (0 deg. az.; 0 deg. el.)
% y-axis is green, point marker denotes poitive y
%   positive y points to (90 deg. az.; 0 deg. el.)
% z-axis is blue, point marker denotes poitive z
%   positive z points to (0 deg. az.; 90 deg. el.)

% Note that you can change the coordinate convention with the 'coord'
% attribute


%% by default your data is interpolated with a very simple algorithm that
%  does not account for the spherical nature of the data, and distorts the
%  data at the north and south pole. If this is to scetchy for you can plot
%  use the 'sph_proc' flag to plot a simple trinangulation that displays
%  the data as it is or use spherical splines for interpolation.
%
%  When triangulating, oddly shaped triangles can be discarded using
%  'triPop'. To change this behaviour set 'triPop' to false (default in AKp)

% we take a subset of our full-spherical data for triangulation and 'triPop'
id = g(:,2)>-45 & g(:,2)<45;

AKf
AKp(l(:,id), 'm1', 'g', g(id,:), 'sph_f', 8000, 'sph_proc', 'tri', 'triPop', 1)

%% this time we use spherical splines of order 1. This takes longer, but
%  looks better. For more information see AKsphSplineInterpDemo.m
  
AKf
AKp(l, 'm1', 'g', g, 'sph_f', 8000, 'sph_proc', 'interpSpline1')

%% there are also spherical plots of TOA, ITD, ILD
h.l = l;
h.r = r;

AKf
AKp(h, 'ild1', 'g', g);

%% moreover, you can also plot spherical data as it is using the 'x' plot
%  type. In this example we plot a single bin of the impulse response just
%  to illustrate how to use this - this plot is pretty useless in
%  particular

AKf
AKp(l(50,:), 'x1', 'g', g, 'cm', 'AKbwr')

%% --------------------------------------------------------- 4. polar plots
%  polarplots are using Matlabs polarplot or mmpolar.m by D.C. Hanselman in
%  earlier matlab versions

% once again we take our horizontal plane HRIRs
id = AKsubGrid(g, 'transverse', 0);

AKf
AKp(l(:,id), 'm6', 'az', g(id,1), 'dr', [-40 20], 'sph_f', 8000)

%% --------------------------------------------- 5. plots of multiple types
%  AKp can also be used to quickly plot multiply things.
%  (see help in AKp and AKpMulti)

% taking a single position
id = AKsubGrid(g, 'any', [90 0]);

% using a predefined plot type
% (you can easily define more plots inside AKpMulti)
AKf
AKp([l(:,id) r(:,id)], '1a')

%% defining a plot yourself
%  in this case the size of the cell array passed as second argument
%  defines the number of subplots - try it yourself:
%  {'tc2d'; 'ms2d'}
%  {'tc2d'  'ms2d'}
%  {'tc2d'  'ms2d'; 'ps2d' []}
AKf
AKp([l(:,id) r(:,id)], {'t2d'; 'm2d'})


%% taking multiple positions

% once again we take our horizontal plane HRIRs
id = AKsubGrid(g, 'transverse', 0);

AKf
AKp(l(:,id), {'t3d'; 'm3d'}, 'y', g(id,1), 'cm', 'AKbwr')


%% ----------------------------------------------------------- 6. HpIR plot
%  there is a quick function that plots
%  - averaged headphone impulse responses (HpIRa)
%  - averaged HpIRs compared to a target function
%  - difference between left and right channel in auditory filters

AKdependencies('FABIAN')

% load HpIRs
HpIRs = SOFAload('HpIRs.sofa');

% load target - we use FABIANs diffuse field HRIR for this
target = SOFAload('FABIAN_CTF_measured_smoothed.sofa');
target = squeeze(target.Data.IR);

AKpHeadphone(HpIRs, target)
