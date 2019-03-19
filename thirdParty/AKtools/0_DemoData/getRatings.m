% This script generates random ratings that mimic results from a listening
% test. They are intended to demonstrate the usage of AKboxplot.m
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
rng(1)

% generate randomly distributed ratings ratings
% 8  dependent variables
% 60 saubjects
% 3  test conditions
ratings = randn(60, 8, 3);
% scale ratings
ratings = AKm(ratings, [8 9 5 12 8 3 12 9], '/');
% change the mean
ratings = AKm(ratings, [.7 -.2 -.5 .5 .8 0 0 .3], '+');
% clip to [-1 1]
ratings = max(ratings, -1);
ratings = min(ratings, 1);