% function [tmp50, tmp95, tmp50_interchannel_mean, tmp95_interchannel_mean, echo_dens] = AKdataBasedMixingTime(IR,N,fs,peak_secure_margin,speedUp)
% Computes the perceptual mixing time from data based predictors
% (echo density) as described in:
%
% Lindau, A.; Kosanke, L.; Weinzierl, S.: 
% "Perceptual evaluation of model- and signal-based predictors of the 
% mixing time in binaural room impulse responses", In:  J. A. Eng. Soc.
%
% see AKperceptualMixingTimeDemo for examples
%
%
% I N P U T:
% IR                 - impulse response of size [N C] where N is thenumber
%                      of samples, and C the number of channels
% N                  - window length in samples (see Abel 2006)
% fs                 - sampling frequency in Hz
% peak_secure_margin - safety margin in samples for onset detection
% speedUp            - stops the time consuming calculation of the echo
%                      density after finding the first value > 1, which is
%                      used for predicting the mixing time. true or false,
%                      default = false.
%
% O U T P U T:
% tmp50                   - average perceptual mixing time
% tmp95                   - 95%-point perceptual mixing time
% tmp50_interchannel_mean - interchannel average of average perceptual mixing time
% tmp95_interchannel_mean - interchannel average of average of 95%-point perceptual mixing time
% echo_dens               - echo density vector calculated inside
%                           AKmixingTimeAbel.m
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
function [tmp50, tmp95, tmp50_interchannel_mean, tmp95_interchannel_mean, echo_dens] = AKdataBasedMixingTime(IR,N,fs,peak_secure_margin,speedUp)

fprintf('\nCalcluating perceptual mixing times tmp50 and tmp95 (in ms) from data-based predictors ...')

% set default value
if ~exist('speedUp', 'var')
    speedUp = false;
end

% preallocate
t_abel      = zeros(1,size(IR,2));
echo_dens   = zeros(length(IR),size(IR,2));

for n = 1:size(IR,2)
    [t_abel(n),echo_dens(:,n)] = AKmixingTimeAbel(IR(:,n),N,fs,peak_secure_margin,speedUp);
end

% tmp from regression equations
tmp50 = 0.8 .* t_abel - 8;
tmp95 = 1.77 .* t_abel -38;


% clip negative values
if sum(tmp50<0)
    idx = find(tmp50<0);
    for g = 1:length(idx)
        tmp50(idx(g)) = 1;
    end
end

if sum(tmp95<0)
    idx = find(tmp95<0);
    for g = 1:length(idx)
        tmp95(idx(g)) = 1;
    end
end

% average perceptual mixing time over all channels
if size(IR,2)>1
    t_abel_interchannel_mean = mean(t_abel);
    tmp50_interchannel_mean = 0.8 .* t_abel_interchannel_mean - 8;
    tmp95_interchannel_mean = 1.77 .* t_abel_interchannel_mean -38;
else
    tmp50_interchannel_mean = [];
    tmp95_interchannel_mean = [];
end

% output the results
fprintf('finished!\n')
if size(IR,2) == 1
    fprintf('tmp50 = %4.2f\n', tmp50)
    fprintf('tmp95 = %4.2f\n', tmp95)
else
    
    fprintf('tmp50 = [')
    for nn = 1:size(IR,2)
        fprintf(' %4.2f', tmp50(nn))
    end
    fprintf(' ], inter channel mean: %4.2f\n', tmp50_interchannel_mean)
      
    fprintf('tmp95 = [')
    for nn = 1:size(IR,2)
        fprintf(' %4.2f', tmp95(nn))
    end
    fprintf(' ], inter channel mean: %4.2f\n', tmp95_interchannel_mean)
end