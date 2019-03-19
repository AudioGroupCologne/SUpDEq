% [N, M, Nsh] = AKgetNM(n)
% returns all combinations of spherical harmonic order n, and degree m up
% to the maximum order N=n. Nsh is the number of spherical harmonics up to
% order N=n.
%
% e.g.
% AKgetNM(1) will output
% N   = 0  1  1  1 (order)
% M   = 0 -1  0 -1 (degree), and
% Nsh = 4
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
function [N, M, Nsh] = AKgetNM(n)

% number of spherical harmonics, N,M,az matrices
Nsh = (n+1)^2;
N   = zeros(Nsh, 1);
M   = N;

% get combinations
for nn = 1:n
    pos_start = (nn)^2+1;   % number of previously calculated n,m-combinations
    N(pos_start:pos_start+2*nn) = nn;
    M(pos_start:pos_start+2*nn) = -nn:nn;
end
