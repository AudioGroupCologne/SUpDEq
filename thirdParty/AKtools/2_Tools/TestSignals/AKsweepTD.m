% [y, N] = AKsweepTD(type, N, f, fs, fade, phaseCorr)
% calculates a sweep in the time domain based on [1], [2]. An optional
% phase correction can be applied for optimal estimatin of non-linearities
% if using the exponential sweep.
%
% See AKtestSignalsDemo.m for examples
%
%
% I N P U T (default)
% type ('exp')      - 'lin' for a sweep with linear increase of frequency
%                     over time
%                     'exp' for a sweep wiht exponential increase of
%                     frequency over time. (Sometimes also called log-sweep
%                     because the group delay increases logarithmically)
% N (2^16)          - sweep length in samples
% f ([50 fs/2])     - start and stop frequency in Hz
% fs (44100)        - sampling frequency in Hz
% fade (20)         - length of fade out in samples. Use 0, or false for
%                     not applying a fade out.
% phaseCorr (false) - correction of the phase according to [1, eq. 16-22].
%                     Note that this changes the length of the sweep.
%
% OUTPUT
% y  - sweep time signal
% N  - length in Samples
%
% [1] Farina, Angelo (2000): "Simultaneous measurement of impulse response
%     and distortion with a swept-sine technique." In: 108th AES Convention
%     , Paris: France
% [2] A. Novak et al. (2010): 'Nonlinear system identification using
%     exponential swept-sine signal.' In: IEEE Trans. Instrumentation and
%     Measurement. 59(8):2220-2229.
%
% v1 03/2013 fabian.brinkmann@tu-berlin.dem Audio Communication Group,
%            TU Berlin

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
function [y, N] = AKsweepTD(type, N, f, fs, fade, phaseCorr)

if ~exist('type', 'var')
    type = 'exp';
end
if ~exist('N', 'var')
    N = 2^16;
end
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('f', 'var')
    f = [50 fs/2];
end
if ~exist('phaseCorr', 'var')
    phaseCorr = false;
end
if ~exist('fade', 'var')
    fade = 20;
end

% time vector and sweep duration
t = (0:N-1)'/fs;
T = t(end);


if strcmpi(type, 'lin')
    
    % get the sweep [1, p. 5]
    y = sin(2*pi*f(1).*t + 2*pi*(f(2)-f(1))/T .* t.^2/2);
    
elseif strcmpi(type, 'exp')
    
    c = log(f(2)/f(1));
    
    if phaseCorr
        
        % f(1)*L has to be an integer for phase correction [2, eq. 22]
        L = 1/f(1) * max(round(T*f(1)/c), 1);
        
        % calculate the T that is actually used, because the rounding operation
        % changes it [2, eq. 14]
        T = L*c;
        N = T*fs;
        t = (0:N-1)'/fs;
        
    else
        
        L = T/c;
        
    end
    
    % get the sweep [1, p. 6], [2, eq. 22]
    y = sin(2*pi*f(1) * L * (exp(t/L)-1));
    
else
    
    error('AKsweepTD:input', 'type must be ''lin'' or ''exp''.')
    
end

% fade out
if fade
    y = AKfade(y, [], 0, fade);
end
