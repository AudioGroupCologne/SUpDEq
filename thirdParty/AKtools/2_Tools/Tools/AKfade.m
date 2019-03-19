% y = AKfade(x, N, Nin, Nout)
%
% applies sine square fade in, and fade out to input signal x.
%
%
% I N P U T:
% x    - Impulse responses [N M C], where N is the number of samples, M the
%        number of measurements, and C the number of channels.
% N    - scalar to set the length of output IRs in samples or two element 
%        vector to specify the start and end point for truncation.
%        Pass [] if you want to keep the
%        original length
% Nin  - length of fade in in samples
% Nout - length of fade out in samples
%
%
% O U T P U T:
% y    - output IRs
%
% v1 04/2016 Fabian Brinkmann, fabian.brinkmann@tu-berlin.de,
%            Audio Communicatin Group, TU Berlin,
%            original dev
% v2 02/2018 Hannes Helmholz, helmholz@campus.tu-berlin.de,
%            Audio Communicatin Group, TU Berlin,
%            fix using also sine (instead of cos) window for fade out
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
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function y = AKfade(x, N, Nin, Nout)

% ------------------------------------------------------------- check input
if isempty(N)
    N = size(x,1);
end

if numel(N)==1
    N = [1 N];
end

N = sort(N);

% ---------------------------------------------------------------- truncate
y = x(N(1):N(2),:,:);

% ----------------------------------------------------------------- fade in
if Nin
    fadeIn       = repmat((sin(linspace(0, pi/2, Nin)).^2)',  [1 size(x,2) size(x,3)]);
    y(1:Nin,:,:) = y(1:Nin,:,:) .* fadeIn;
end

% ---------------------------------------------------------------- fade out
if Nout
    fadeOut               = repmat((sin(linspace(pi/2, 0, Nout)).^2)', [1 size(x,2) size(x,3)]);
    y(end-Nout+1:end,:,:) = y(end-Nout+1:end,:,:) .* fadeOut;
end
