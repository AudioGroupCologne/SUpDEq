% y = AKedc(x, norm)
%
% calculates normalized (norm = 1, default) or unnormalized (norm = 0)
% energy decay curve (EDC) in dB of input data x [samples x channels]
%
% See AKplotDemo.m and AKp.m for examples
%
% F. Brinkmann, Audio Communicatin Group TU-Berlin
% fabian.brinkmann@tu-berlin.de, 04/2012

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
function y = AKedc(x, norm)

if ~exist('norm', 'var')
    norm = 1;
end

y = flipud(cumtrapz(flipud(x.^2)));

if norm
    y = y./repmat(trapz(x.^2), size(x,1),1);
end

y = 10*log10(y);
