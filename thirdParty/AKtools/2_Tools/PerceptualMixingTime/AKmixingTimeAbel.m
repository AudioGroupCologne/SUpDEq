% [t_abel,echo_dens] = AKmixingTimeAbel(IR,N,fs,peak_secure_margin,speedUp)
% Computes the transition time between early reflections and stochastic
% reverberation based on the assumption that sound pressure in a 
% reverberant field is Gaussian distributed. From: 
%
% Abel & Huang 2006, "A simple, robust measure of reverberation echo
% density", In: Proc. of the 121st AES Convention, San Francisco
%
% see AKperceptualMixingTimeDemo for examples
%
%
% I N P U T:
% IR                 - impulse response (1 channel only!)
% N                  - window length in Samples
% fs                 - sampling rate in Hz
% peak_secure_margin - safety margin in samples for onset detection
% speedUp            - stops the time consuming calculation of the echo
%                      density after finding the first value > 1, which is
%                      used for predicting the mixing time. true or false,
%                      default = false.
%
% O U T P U T:
% t_abel             - mixing time after Abel & Huang (2006, echo density = 1)
% echo_dens          - echo density vector
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
function [t_abel,echo_dens] = AKmixingTimeAbel(IR,N,fs,peak_secure_margin,speedUp)

% set default value
if ~exist('speedUp', 'var')
    speedUp = false;
end

% preallocate
s           = zeros(1,length(IR));
anz_s       = zeros(1,length(IR));
echo_dens   = nan(1,length(IR));

if length(IR) < N
    error('IR shorter than analysis window length (1024 samples). Provide at least an IR of some 100 msec.')
end

% value for normalizing the mixing time
echo_dens_norm = 1 / erfc(1/sqrt(2));

for n = 1:length(IR)
    % window at the beginning (increasing window length)
    if n <= N/2+1
        
        % standard deviation 
        s(n)        = std(IR(1:n+N/2-1));
        
        % number of tips outside the standard deviation
        anz_s(n)    = sum(abs(IR(1:n+N/2-1))>s(n)); 
        
        % echo density
        echo_dens(n)= anz_s(n)/N;                                   
    end

    % window in the middle (constant window length)
    if n > N/2+1 && n <= length(IR)-N/2+1    
        s(n)        = std(IR(n-N/2:n+N/2-1));
        anz_s(n)    = sum(abs(IR(n-N/2:n+N/2-1))>s(n));
        echo_dens(n)= anz_s(n)/N;                            
    end

    % window at the end (decreasing window length)
    if n > length(IR)-N/2+1   
        s(n)        = std(IR(n-N/2:length(IR)));
        anz_s(n)    = sum(abs(IR(n-N/2:length(IR)))>s(n));
        echo_dens(n)= anz_s(n)/N;
    end
    
    % normalize echo density
    echo_dens(n) = echo_dens(n) * echo_dens_norm;       
        
    if echo_dens(n) > 1 && speedUp
        break
    end
    
end                   

% transition point (Abel & Huang (2006))
% (echo density first time greater than 1)
d           = min(find(echo_dens>1));    %#ok<MXFND>
t_abel      = (d-peak_secure_margin)/fs*1000;

if isempty(t_abel)
    fprintf('\n')
    error('Mixing time not found within given temporal limits. Try again with extended stopping crtiterion.')
end