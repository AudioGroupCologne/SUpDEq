% [e_n, N] = AKch(n, phi, mode)
% calculates circular harmonics (CH) basis functions according to [1],
% equation (4).
%
% Example usage:
% See AKcircularHarmonicsDemo.m for examples
%
% INPUT
% n    - CH order. All CH for order -n to b will be calcualated
% phi  - angular positions in degree. Vector of size [Q x 1] with values
%        between 0 and 260 degree.
% mode - 'complex' or 'real' for calculating complex (default) and real
%        valued CH. In the latter case the imaginary part of the complex 
%        valued CH is taken for n<0 and the real part for n>0.
%
% OUTPUT
% e_n - circular harmonics of size [Q x (2*N)+1]. ALL CH for one phi-value
%       are saved in a row. Columns specifiy order n.
% N   - CH Order corresponding to each column of e_n
%
% [1] Noam R Shabtai, Gottfried Behler, and Michael Vorlaender: "Generation
%     of a reference radiation pattern of string instruments using
%     automatic excitation and acoustic centering." J. Acoust. Soc. Am.,
%     138(5): EL479-EL456, (2015).
%
% fabian.brinkmann@tu-berlin.de,
% Audio Communication Group, TU Berlin
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
function [e_n, N] = AKch(n, phi, mode)

% --- input parsing ---
if numel(n)~=1
    error('AKsh:Input', 'n must be a scalar value')
elseif rem(n,1)
    error('AKsh:Input', 'n must be an integer value')
end

if~exist('mode', 'var')
    mode = 'complex';
end

phi = reshape(phi, [numel(phi) 1]);

% Combnation of angles and orders
phi  = phi/180*pi;
N   = -n:n;
NN  = AKm(N, phi, '*');

% calculate circular harmonics
e_n = exp(1j*NN) / (2*pi);

% get real valued circular harmonics
if strcmpi(mode, 'real')
    e_n(:, N< 0) = imag(e_n(:, N< 0));
    e_n(:, N>=0) = real(e_n(:, N>=0));
end
