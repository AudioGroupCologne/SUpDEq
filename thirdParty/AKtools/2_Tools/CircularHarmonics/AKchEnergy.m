% [e_n, avg_n, e_N, avg_N] = AKchEnergy(d_n)
% calcualtes the energy of the circular function f(phi) using the
% corresponding circular harmonics coefficients d_n. The energy is
% calculated using Parseval's theorem according to Eq. (1.43) in [1] and
% normalized by 1/(4*pi^2).
%
% I N P U T
% d_n    - cicrular harmonics coefficients of size [(N+1)*2 x M], where N
%          is the circular harmonics order and M the number of sets of
%          coefficients, e.g., frequency bins.
%
% O U T P U T
% e_n    - energy for each set of d_n. Size [1 x M]
% avg_n  - energetic average of e_n across M. Scalar value.
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
function [e_n, avg_n, e_N, avg_N] = AKchEnergy(d_n)

% spherical harmonics order
N = ( size( d_n, 1 ) - 1 ) / 2;
F = size(d_n, 2);

% energy contained in f_nm coefficients for each bin (diffuse field transfer function)
e_n = 1/(4*pi^2) * sum( abs(d_n).^2, 1 );

% energetic average across all bins
avg_n = sum( e_n, 2 ) / F;

% energey contained within orders
e_N = zeros(N+1, F);

for nn = 0:N
    
    id = unique( [N+1+nn N+1-nn] );
    
    % energy contained in current order
    e_N(nn+1,:,:) = 1/(4*pi^2) * sum( abs( d_n(id,:,:) ).^2, 1 );
end

% energetic average within orders and across bins
avg_N = sum( e_N, 2 ) / (N+1);