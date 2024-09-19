% y = AKicht(d_n, doIFFT, phi, CHTmode, is_even, compact, CHmode)
% does a discrete inverse spherical harmonics/fourier transform according
% to eq. (3) in [1].
%
% See AKcircularHarmonicsDemo.m for examples
%
% I N P U T:
% d_n     - circular harmonics coefficients with formar specified by AKcht
% doIFFT  - if true (default) an inverse Fourier transform is done after
%           the inverse circular harmonics transform
% phi     - the spatial sampling grid in form of a Qx1 vector holding
%           the angular values. Coordinate convention according to AKch.m
% CHTmode - as specified by AKcht (default='db_unwrap')
% is_even - as specified by AKcht (default=1)
% compact - as specified by AKcht (auto detected)
% CHmode  - as specified by AKcht (default = 'complex')
%
% O U T P U T:
% y       - reconstructed spatial function at points specified by phi.
%           Each column holds the data for one point of the spatial grid.
%
% [1] Noam R Shabtai, Gottfried Behler, and Michael Vorlaender: "Generation
%     of a reference radiation pattern of string instruments using
%     automatic excitation and acoustic centering." J. Acoust. Soc. Am.,
%     138(5): EL479-EL456, (2015).
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
function y = AKicht(d_n, doIFFT, phi, CHTmode, is_even, compact, CHmode)

if ~exist('CHTmode', 'var')
    CHTmode = 'db_unwrap';
end
if ~exist('is_even', 'var')
    is_even = 1;
end
if ~exist('compact', 'var')
    if size(d_n,3) == 1 && any(strcmpi(CHTmode, {'db_unwrap' 'abs_unwrap' 'real_imag'}))
        compact = true;
    else
        compact = false;
    end
end
if ~exist('CHmode', 'var')
    CHmode = 'complex';
end

phi = reshape(phi, [numel(phi), 1]);

% ---------------------------- get e_m matrix if not passed to the function
N   = ( size(d_n,1) -1 ) / 2;
e_n = AKch(N, phi(:,1), CHmode);

% ---------------------------------------------------- inverse CH transform
Y = zeros(size(e_n,1), size(d_n,2), size(d_n,3));

for n = 1:size(d_n,2)
    for c = 1:size(d_n,3)
        Y(:,n,c) = e_n * d_n(:,n,c);
    end
end

% transpose (could have been done above, but this makes it more
% illustrative)
tmp = Y;
Y   = zeros(size(d_n,2), size(e_n,1), size(d_n,3));
for c = 1:size(d_n,3)
    Y(:,:,c) = tmp(:,:,c).';
end

clear tmp

% -------------------------------- post process x to match desired CHT mode
if compact && any(strcmpi(CHTmode, {'db_unwrap' 'abs_unwrap' 'real_imag'}))
    Y   = cat(3, real(Y), imag(Y));
end

switch CHTmode
    case 'db_unwrap'
        y = 10.^(Y(:,:,1)/20) .* exp(1j*Y(:,:,2));
    case 'abs_unwrap'
        y = Y(:,:,1) .* exp(1j*Y(:,:,2));
    case 'complex'
        y = Y;
    case 'real_imag'
        y = Y(:,:,1) + 1j*Y(:,:,2);
    case 'db'
        y = 10.^(Y/20);
    case 'abs'
        y = Y;
    otherwise
        error(['CHTmode ''' CHTmode ''' unknown'])
end

% ------------------------------------ transfer into time domain if desired
if doIFFT
    y = ifft(AKsingle2bothSidedSpectrum(y, is_even), 'symmetric');
end
