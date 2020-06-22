% y = AKisht(f_nm, doIFFT, multi, SHTmode, is_even, compact, SHmode)
% does a discrete inverse spherical harmonics/fourier transform according
% to eq. (3.34) in [1].
%
% See AKsphericalHarmonicsTransformDemo.m for examples
%
% I N P U T:
% f_nm    - spherical harmonics coefficients with formar specified by AKsht
% doIFFT  - if true (default) an inverse Fourier transform is done after
%           the inverse spherical harmonics transform
% multi   - either the desired spatial coordinates (see help AKsh for
%           coordinate convention) or the Ynm matrix holding the values of
%           the spherical harmonics functions. The sampling grid must be
%           passed as a Qx2 matrix where the first column holds the azimuth
%           and the second the elevation.
% SHTmode - as specified by AKsht (default='db_unwrap')
% is_even - as specified by AKsht (default=1)
% compact - as specified by AKsht (auto detected)
% SHmode  - as specified by AKsht (default = 'complex')
%
% O U T P U T:
% y       - reconstructed spatial function at points specified by multi.
%           Each column holds the data for one point of the spatial grid.
%
% [1] Boaz Rafaely: Fundamentals of spherical array processing. In.
%     Springer topics in signal processing. Benesty, J.; Kellermann, W.
%     (Eds.), Springer, Heidelberg et al. (2015).
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
% 09/2015

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
function y = AKisht(f_nm, doIFFT, multi, SHTmode, is_even, compact, SHmode)

if ~exist('SHTmode', 'var')
    SHTmode = 'db_unwrap';
end
if ~exist('is_even', 'var')
    is_even = 1;
end
if ~exist('compact', 'var')
    if size(f_nm,3) == 1 && any(strcmpi(SHTmode, {'db_unwrap' 'abs_unwrap' 'real_imag'}))
        compact = true;
    else
        compact = false;
    end
end
if ~exist('SHmode', 'var')
    SHmode = 'complex';
end

% ---------------------------- get Ynm matrix if not passed to the function
if size(multi, 2)==2
    N   = sqrt(size(f_nm,1))-1;
    Ynm = AKsh(N, [], multi(:,1), multi(:,2), SHmode);
else
    Ynm = multi;
end

% ---------------------------------------------------- inverse SH transform
Y = zeros(size(Ynm,1), size(f_nm,2), size(f_nm,3));

for n = 1:size(f_nm,2)
    for c = 1:size(f_nm,3)
        Y(:,n,c) = Ynm * f_nm(:,n,c);
    end
end

% transpose (could have been done above, but this makes it more
% illustrative)
tmp = Y;
Y   = zeros(size(f_nm,2), size(Ynm,1), size(f_nm,3));
for c = 1:size(f_nm,3)
    Y(:,:,c) = tmp(:,:,c).';
end

clear tmp

% -------------------------------- post process x to match desired SHT mode
if compact && any(strcmpi(SHTmode, {'db_unwrap' 'abs_unwrap' 'real_imag'}))
    Y   = cat(3, real(Y), imag(Y));
end

switch SHTmode
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
        error(['SHTmode ''' SHTmode ''' unknown'])
end

% ------------------------------------ transfer into time domain if desired
if doIFFT
    y = ifft(AKsingle2bothSidedSpectrum(y, is_even), 'symmetric');
end
