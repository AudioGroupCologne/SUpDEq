% [e_nm, avg_nm, e_N, avg_N] = AKshEnergy(f_nm)
% calcualtes the energy of the sphercial function f(phi, theta) using the
% corresponding spherical harmonics coefficients f_nm. The energy is
% calculated using Parseval's theorem according to Eq. (1.43) in [1] and
% normalized by 1/(4*pi).
%
% I N P U T
% f_nm   - spherical harmonics coefficients of size [(N+1)^2 x M], where N
%          is the spherical harmonics order and M the number of sets of
%          coefficients, e.g., frequency bins.
%
% O U T P U T
% e_nm   - energy for each set of f_nm. Size [1 x M]
% avg_nm - energetic average of e_nm across M. Scalar value.
% e_N    - energy contained in each order. Size [N x M]
% avg_N  - energetic average of e_N across N. Size [N x 1]
%
%
% [1] Boaz Rafaely: Fundamentals of spherical array processing. In.
%     Springer topics in signal processing. Benesty, J.; Kellermann, W.
%     (Eds.), Springer, Heidelberg et al. (2015).
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
% 02/2018

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
function [e_nm, avg_nm, e_N, avg_N] = AKshEnergy(f_nm)

% spherical harmonics order
N = sqrt( size( f_nm, 1 ) ) - 1;
F = size(f_nm, 2);

% energy contained in f_nm coefficients for each bin (diffuse field transfer function)
e_nm = 1/(4*pi) * sum( abs(f_nm).^2, 1 );

% energetic average across all bins
avg_nm = sum( e_nm, 2 ) / F;

% energey contained within orders
e_N = zeros(N+1, F);

for nn = 0:N
    pos_start = (nn)^2+1;                   % number of previously averaged coefficients
    id        = pos_start:pos_start+2*nn;   % indieces of coefficents belonging to the current order
    
    % average the current order
    e_N(nn+1,:,:) = 1/(4*pi) * sum( abs( f_nm(id,:,:) ).^2, 1 );
end

% energetic average within orders and across bins
avg_N = sum( e_N, 2 ) / (N+1);