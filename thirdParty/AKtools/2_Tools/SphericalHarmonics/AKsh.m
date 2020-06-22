% [Ynm, N, M] = AKsh(n, m, az, el, mode)
% calculates complex spherical harmonics (SH) basis functions according
% to [1], equations (1.9, 1.30, 1.34). Real valued SH functions are
% caluclated according to [2].
%
% Example usage:
% [Ynm, N, M] = AKsh(3, [], [0 90 180 270], [90 90 90 90]);
%   computes all SH up to order 3 at for positions on the horizontal plane
% [Ynm, N, M] = AKsh(5, -2, [0 90 180 270], [90 90 90 90]);
%   computes SH of order 3 and degree -2 at same positions
%
% See AKsphericalHarmonicsDemo.m for examples
%
% INPUT
% n    - SH order. All SH up to order n will be calcualated if m is empty,
%        only SH according to n, and m, otherwise
% m    - SH degree, see above.
% az   - azimuth positions in degree. Can be scalar or vector of size(el)
%        (0=front, 90=left, 180=back, 270=right)
%        (0 points to positive x-axis, 90 to positive y-axis)
% el   - elevations in degree. Can be scalar or vector of size(az)
%        (0=North Pole, 90=front, 180= South Pole)
%        (0 points to positive z-axis, 180 to negative z-axis)
% mode - 'complex' or 'real' for calculating complex (default) and real
%        valued SH. In the latter case the imaginary part of the complex 
%        valued SH is taken for m<0 and the real part for m>0.
%
% OUTPUT
% Ynm - Complex valued spherical harmonics. One SH for all (az,el) is saved
%       in each Column. Consequently columns specifiy order n, and degree m
%       and rows specify (az, el).
% N   - SH Order corresponding to Ynm
% M   - SH Degree corresponding to Ynm
%
%
% [1] Boaz Rafaely: Fundamentals of spherical array processing. In.
%     Springer topics in signal processing. Benesty, J.; Kellermann, W.
%     (Eds.), Springer, Heidelberg et al. (2015).
% [2] Earl G. Williams: Fourier Acoustics. Sound radiation and nearfield
%     acoustical holography.  Academic Press, San Diego et al., (1999).
% [3] Franz Zotter: Analysis and synthesis of sound-radiation with 
%     spherical arrays. Ph.D. dissertation, University of Music and
%     Performing arts (2009).
%
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
function [Ynm, N, M] = AKsh(n, m, az, el, mode)

% --- input parsing ---
if numel(n)~=1
    error('AKsh:Input', 'n must be a scalar value')
elseif rem(n,1)
    error('AKsh:Input', 'n must be an integer value')
end

if~exist('mode', 'var')
    mode = 'complex';
end

az = reshape(az, [numel(az) 1]);
el = reshape(el, [numel(el) 1]);

if size(az,1)~=size(el,1) || size(az,2)~=size(el,2)
    error('AKsh:Input', 'az and el must be of same size')
end

% Number of angles and coordinate conversion
Naz = numel(az);
az  = az/180*pi;
el  = el/180*pi;

% --- get combinations of n and m to compute ---
if isempty(m) || ~isnumeric(m)  % --> we compute all m for nn = 0;n
    
    % number of spherical harmonics, N,M,az matrices    
    [N, M, Nsh] = AKgetNM(n);
    
    N  = repmat(N', [numel(az) 1]);
    M  = repmat(M', [numel(az) 1]);
    az = repmat(az, [1 Nsh]);
    
else   % --> we compute one n,m-combination
    % input check
    if numel(m)~=1
        error('AKsh:Input', 'm must be a scalar value')
    elseif rem(m,1)
        error('AKsh:Input', 'm must be an integer value')
    elseif abs(m)>n
        error('AKsh:Input', '|m| must not be larger than n')
    end
    
    % N,M-matrices
    Nsh = 1;
    N   = repmat(n, [Naz Nsh]);
    M   = repmat(m, [Naz Nsh]);
end

clear nn num_sh pos_start Naz

% --- calculate spherical harmonics ---
% scalar (dependend on n,m; independend from az, el)

% check if all factorials can be computed
if any( [N-M N+M] > 170 )
    error('AKist:input', 'SH basis functions can only be computed up to order 85. Otherwise the factorials can not easiliy be computed with built in variable types.')
end

if strcmpi(mode, 'complex')
    a = sqrt((2*N+1)./(4*pi) .* factorial(N-M)     ./factorial(N+M));
else
    a = sqrt((2*N+1)./(4*pi) .* factorial(N-abs(M))./factorial(N+abs(M)));
end
% azimuthal changing part of spehrical harmonics
if strcmpi(mode, 'complex')
    e = exp(1j*M.*az);
else
    e = ones(size(M));
    id    = M<0;
    e(id) = sqrt(2) * sin(M(id).*az(id));
    id    = M>0;
    e(id) = sqrt(2) * cos(M(id).*az(id));
end
% elevation dependend part of spherical harmonics
if Nsh > 1
    l = ones(size(N));
    % get combinations
    for nn = 1:n
        pos_start = (nn)^2+1;   % number of previously calculated n,m-combinations
        
        % calculate legendre functions for positive m
        l(:,pos_start+nn:pos_start+2*nn) = legendre(nn, cos(el))';
        if strcmpi(mode, 'complex')
            % legendre functions for negative m from [2], eq. (6.31)
            mm = 1:nn;
            lScale = (-1).^mm .* factorial(nn-mm)./factorial(nn+mm);
            l(:,pos_start:pos_start+nn-1) = fliplr(l(:,pos_start+nn+1:pos_start+2*nn) .* repmat(lScale, [size(az,1) 1]));
        else
            l(:,pos_start:pos_start+nn-1) = fliplr( l(:,pos_start+nn+1:pos_start+2*nn) );
        end
    end
else
    l = legendre(n, cos(el))';
    if ~isempty(m)
        l = l(:,abs(m)+1);
        if m < 0 && strcmpi(mode, 'complex')
            % legendre function for negative m from [2], eq. (6.31)
            l = l * (-1).^-m .* factorial(n+m)./factorial(n-m);
        end
    end
end

if strcmpi(mode, 'real')
    % remove Condon-Shortley phase
    l = (-1).^M .* l;
end

% get spherical harmonics
Ynm = a .* l .* e;

% get corresponding n,m-vectors
N = N(1,:);
M = M(1,:);
