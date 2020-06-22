% sine = AKsine(T, f, fs, zeroEnd)
% 
% generates sine signal(s)
% See AKtestSignalsDemo.m for examples
%
% I N P U T (default value):
% T (5)                - duration of sine signal in [s]
% f (1000)             - frequency [Hz]. Scalar or vector. If f is a vector
%                        sine will have as many channels as f has elements.
% fs (44100)           - sampling frequency [Hz]
% zeroEnd (true)       - true: T is modified so that sine ends with a zero,
%                              i.e. with a phase of 0 or 180 degree.
%                        false: sine might end with non-zero element.
%
% v1 12/2015 fabian.brinkmann@tu-berlin.de, Audio Communication Group,
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
function sine = AKsine(T, f, fs, zeroEnd)

% set default parameter
if ~exist('T', 'var')
    T = 5;
end
if ~exist('f', 'var')
    f = 1000;
end
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('zeroEnd', 'var')
    zeroEnd = true;
end

% initital duration in samples
Ninit = ceil(T*fs);

% adjust duration
if zeroEnd
    
    % period duration of  frequencies
    Nf = 1./f * fs;
    
    % missing samples to end at 0, or 180 degree phase shift
    dN = Nf - rem(Ninit, Nf);
    
    % final duration in samples (append one zero in case we have a period
    % length that is not of an integer sample duration)
    N = Ninit + ceil(max(dN))+1;
    
    % readjust missing samples to final sine duration
    dN = round(N - rem(N, Nf));
    
else
    
    N  = Ninit;
    dN = N * ones(numel(f), 1);
    
end

% allocate space
sine = zeros(N, numel(f));

% time vector
t    = (0:1/fs:1/fs*(N-1))';

% get sine signals
for cc = 1:numel(f)
    sine(1:dN(cc),cc) = sin(2*pi*f(cc)*t(1:dN(cc)));
end
