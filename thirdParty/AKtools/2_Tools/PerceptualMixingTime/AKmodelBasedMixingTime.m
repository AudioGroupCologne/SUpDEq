% [tmp50, tmp95] = AKmodelBasedMixingTime(h,b,l)
% Computes the perceptual mixing time from model based predictors as
% described in:
%
% Lindau, A.; Kosanke, L.; Weinzierl, S.: 
% "Perceptual evaluation of model- and signal-based predictors of the 
% mixing time in binaural room impulse responses", In:  J. A. Eng. Soc.
%
% see AKperceptualMixingTimeDemo for examples
% 
% I N P U T:
% tmp50    - average perceptual mixing time
% tmp95    - 95%-point perceptual mixing time
%
% O U T P U T:
% h        - height [m]
% b        - width [m]
% l        - length [m]
%
%
% A. Lindau, L. Kosanke, 2011
% alexander.lindau@tu-berlin.de
% audio communication group
% Technical University of Berlin

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
%-------------------------------------------------------------------------%
function [tmp50, tmp95] = AKmodelBasedMixingTime(h,b,l)

% calculate room properties
V       = h * l * b;                % volume
S       = 2*l*b + 2*b*h + 2*l*h;    % surface area
    
fprintf('\nVolume: %5.2f cubic meters.',V)
fprintf('\nSurface area: %5.2f square meters.\n',S)

% predict tmp from linear models 
tmp50  = 20.08 * V/S + 12;
tmp95  = 0.0117 * V + 50.1;

fprintf('\nPerceptual mixing times in ms from model-based predictors:\n')
fprintf('tmp50 = %4.2f\n', tmp50)
fprintf('tmp95 = %4.2f\n', tmp95)


