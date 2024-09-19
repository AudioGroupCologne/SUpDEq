% Demonstration spherical spline interpolation. This might be usefull in
% case you have spherical data that can not be interpolated using spherical
% harmoics.
%
% 11/2014 - Tobias Muenzer, txmunzerx@gmail.com
%
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

% generate spatial sampling grid
g = AKgreatCircleGrid(-90:30:90, 30, 90);

% generate random data
rng(1)
data = rand(size(g,1), 1);

% generate spatial sampling grid for interpolation
g_ref = AKgreatCircleGrid(-90:5:90, 5, 90);

% set interpolation parameters
m          = [1 2 3];       % spline order m={1, 2, 3}
lambda     = [0 .01 0.02];  % smoothing factor
do_plot    = 1;             % 0: no plot, 1: planar plot, 2: spherical plot

%% interpolate with different spline orders and smoothing factors
%  (the colored dots in the plots show the interpolated data, the crosses
%  the position of the original data.)
AKf
for mm = 1:numel(m)
    for ll = 1:numel(lambda)
        
        subtightplot(3,3, (m(mm)-1)*3+ll, [0 0])
        
        data_interp = AKsphSplineInterp(g(:,1), g(:,2), data, g_ref(:,1), g_ref(:,2), m(mm), lambda(ll), 'deg', do_plot);
        
        set(gca, 'xTick', 0:90:270, 'yTick', -45:45:45)
        title(''); box on; grid on
        if ll ~= 1; ylabel ''; set(gca, 'yTickLabel', []); end
        if mm ~= 3 ; xlabel ''; set(gca, 'xTickLabel', []); end
        
        text(5, 80, ['m=' num2str(m(mm)) ', lambda=' num2str(lambda(ll))], 'backgroundColor', 'w')
        
    end
end
