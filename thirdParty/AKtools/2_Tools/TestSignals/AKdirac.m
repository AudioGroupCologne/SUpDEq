% dirac = AKdirac(samples, channels, delay)
% creates dirac signal(s)
%
% See AKtestSignalsDemo.m for examples
%
%
% I N P U T
% samples  - length of dirac signal in samples (default = 1024)
% channels - number of channels (default = 1)
% delay    - delay in samples, scalar or vector with channels elements
%            (default = 0)
%
% O U T P U T
% dirac    - dirac signals of size [sampples x channels]
%
% v1 2015/03 fabian.brinkmann@tu-berlin.de, Audio Communicaiton Group,
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
function dirac = AKdirac(samples, channels, delay)

if ~exist('samples', 'var')
    samples = 1024;
end

if ~exist('channels', 'var')
    channels = 1;
end

if ~exist('delay', 'var')
    delay = 0;
end
if numel(delay) ~= 1 && numel(delay)~= channels
    error('AKsignalDirac:Input', 'delay must be scalar or match the number of channels')
end

dirac = zeros(samples, channels);

if numel(delay) == 1
    dirac(delay+1,:) = 1;
else
    for c=1:channels
        dirac(delay(c)+1,c) = 1;
    end
end
