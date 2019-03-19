% [gnm, D] = AKshRotate(fnm, rot)
% rotates function given by spherical harmonics (SH) coefficients fnm using
% rot in degree. SH convention and rotation follows AKsh, and [1].
%
% See AKsphericalHarmonicsDemo.m for examples
%
% I N P U T:
% fnm - spherical harmonics coefficients in the format specified by AKsht.
% rot - three elements vector specifying successive counter clockwise
%       rotations for rot(3) degree about z-axis, rot(2) degree about
%       y-axis, and rot(1) degree about z-axis
%
% O U T P U T:
% gnm - rotated coefficients fnm
% D   - rotation Matrix gnm(:,i,j)=D*fnm(:,i,j)
%
%
% [1] Boaz Rafaely: Fundamentals of spherical array processing. In.
% Springer topics in signal processing. Benesty, J.; Kellermann, W. (Eds.),
% Springer, Heidelberg et al., first edition (2015).
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
function [gnm, D] = AKshRotate(fnm, rot)

% get rotation angles in radians
rot   = rot/180*pi;
alpha = rot(1);
beta  = rot(2);
gamma = rot(3);

N   = size(fnm,1);  % number of coefficients
F   = size(fnm,2);  % number of columns (e.g. freuquencies)
C   = size(fnm,3);  % number of pages (e.g. coefficient sets or magnitude and phase)

Nsh = sqrt(N)-1;    % SH order

% allocate memory
D = zeros(N);
gnm = zeros(size(fnm));

for n = 0:Nsh
    
    % all degrees for current order
    mm = (-n:n)';
    
    for m = -n:n
        
        % indices of current degrees m, and order n
        id_m  = AKnm2id(n, m);
        % indices of all degrees m for current order n
        id_mm = AKnm2id(n, mm);
        
        % constants mu, nu, s, zeta, eq. (1.74-1.75)
        mu = abs(m - mm);
        nu = abs(m + mm);
        s  = n - (mu+nu)/2;
        zeta       = ones(2*n+1,1);
        zeta(mm<m) = (-1).^(mm(mm<m)-m);
        
        % get jacobi polynomial values used in eq. (1.74)
        jac = zeros(size(mm));
        for j = 1:numel(jac)
            tmp = j_polynomial(1, s(j), mu(j), nu(j), cos(beta));
            jac(j) = tmp(s(j)+1);
        end
        
        % get Wigner-d function values from eq. (1.74)
        Wigner_d = zeta .* sqrt(factorial(s).*factorial(s+mu+nu)./(factorial(s+mu).*factorial(s+nu))) .* ...
            sin(beta/2).^mu .* cos(beta/2).^nu .* jac;
        
        % get Wigner-D function values, eq. (1.73)
        Wigner_D = exp(-1j*m*alpha) * Wigner_d .* exp(-1j*mm*gamma);
        
        % apply rotation, eq. (1.77)
        % (is not done this way but at the end using the D matrix)
        % gnm(id_m) = sum(fnm(id_mm) .* Wigner_D);
        
        % write current Wigner D values in rotation matrix, eq. (1.78)
        D(id_m, id_mm) = Wigner_D;
    end
end

% apply rotation, eq. (1.78)
for c = 1:C
    for f = 1:F
        gnm(:,f,c) = D * fnm(:,f,c);
    end
end
