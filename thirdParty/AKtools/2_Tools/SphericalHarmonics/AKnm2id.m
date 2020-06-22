% id = AKnm2id(n,m)
% returns index id of spherical harmonics (SH) base function of order n,
% and degree m. SH functions can be caclulated with AKsh.m
%
% See AKsphericalHarmonicsDemo.m for an example
%
% I N P U T:
% n, and m give the order and degree of the SH base functions, and must be
% vectors of same size, or one scalar value and one vector.
%
% fabian.brinkmann@tu-berlin.de,
% Audio Communication Group, TU Berlin
% 04/2015

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
function id = AKnm2id(n,m)

n = reshape(n, [numel(n) 1]);
m = reshape(m, [numel(m) 1]);

if numel(n) < numel(m)
    n(2:numel(m),1) = n(1);
elseif numel(m) < numel(n)
    m(2:numel(n),1) = m(1);
end

id = n.^2 + n+1 + m;
