% Demonstrates basic usiage of the spherical harmonics processing. Function
% definition follows [1].
%
% 1. plot all SH basis function up to order N
% 2. plot single SH basis function
% 3. rotation in the SH domain
%
% [1] Boaz Rafaely: Fundamentals of spherical array processing. In.
% Springer topics in signal processing. Benesty, J.; Kellermann, W. (Eds.),
% Springer, Heidelberg et al., first edition (2015).
%
% 04/2015   - fabian.brinkmann@tu-berlin.de,

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
close all; clear; clc

%% -------------------------- 1. plot all spherical harmonics up to order n
%  similar to Fig. 1.5-1.8 in [1]
clear; clc

% spatial sampling grid for plotting
g      = AKgreatCircleGrid(90:-7.5:-90, 7.5, 90, 0);
% definition of elevation is different for spherical harmonics:
% 0 degree: north pole, 90 degree: front, 180 degree south pole
g(:,2) = 90-g(:,2);

% set order
N = 1;
M = [];

% get spherical harmonics basis functions
[Ynm, n, m] = AKsh(N, M, g(:,1), g(:,2), 'real');

% view onto specified axis {'x', 'y', 'z', 'xyz'}
viewPoint = 'xyz';

% plot
AKf(10*(2*N+1),10*N)
for k = 1:numel(n)
    subtightplot(N+1,2*N+1, N+1, [0 0])                         % <- nice function for subplotting!
    subtightplot(N+1,2*N+1, n(k)*(2*N+1) + N+1+ m(k), [0 0])

    if m(k)<0
        % plot imaginary part for m<0
        AKp(imag(Ynm(:,k))', 'x7', 'az', g(:,1), 'el', g(:,2), 'cb', 0, 'sph_proc', 'tri', 'coord', 2)
    else
        % plot real part for m>=0
        AKp(real(Ynm(:,k))', 'x7', 'az', g(:,1), 'el', g(:,2), 'cb', 0, 'sph_proc', 'tri', 'coord', 2)
    end
    
    switch viewPoint
        case 'x'; view([90 0])
        case 'y'; view([180 0])
        case 'z'; view([0 90])
        case 'xyz'
    end
    
    title('')
    pause(-.01)
end

%% ------------------------------------ 2. plot specific spherical harmonic
%  similar to Fig. 1.4 in [1]
clear; clc

% grid for plotting
g      = AKgreatCircleGrid(90:-2:-90, 2, 90, 0);
g(:,2) = abs(g(:,2)-90);

% order and degree
N = 5;
M = -5;

% get spherical harmonic function
Ynm = AKsh(N, M, g(:,1), g(:,2));

% plot it
AKf(40,16)
subtightplot(1,3,1, [0 0])
subtightplot(1,3,1, [0 0])
    AKp(imag(Ynm)', 'x7', 'az', g(:,1), 'el', g(:,2), 'cb', 0, 'sph_proc', 'tri', 'coord', 2)
subtightplot(1,3,2, [0 0])
    AKp(abs(Ynm)' , 'x7', 'az', g(:,1), 'el', g(:,2), 'cb', 0, 'sph_proc', 'tri', 'coord', 2)
subtightplot(1,3,3, [0 0])
    AKp(real(Ynm)', 'x7', 'az', g(:,1), 'el', g(:,2), 'cb', 0, 'sph_proc', 'tri', 'coord', 2)

%% ---------------------- 3. Rotation of a truncated spherical cap function
%  similar to Fig. 1.15 in [1]
clear; clc

% grid for plotting
g      = AKgreatCircleGrid(90:-5:-90, 5, 90, 0);
g(:,2) = abs(g(:,2)-90);

% maximum order
N = 2;

% rotation angles according to [1], eq. (1.77)
rotation_angles = [0 45 0];

% get SH functions
[Ynm, n, m] = AKsh(N, [], g(:,1), g(:,2));

% get f_nm weights [1], eq. (1.63)
fnm = zeros(numel(m), 1);

SHfunc = 'cap';
switch SHfunc
    case 'cap'
        % size of spherical cap [1], p. 22
        alpha = 30/180*pi;
        
        fnm(1) = sqrt(pi)*(1-cos(alpha));
        for nn = 1:N
            Pm1 = legendre(nn-1, cos(alpha));
            Pp1 = legendre(nn+1, cos(alpha));
            fnm(AKnm2id(nn,0)) = sqrt(pi/(2*nn+1))*(Pm1(1)-Pp1(1));
        end
        
        clear nn Pm1 Pp1
    case 'dirac'
        % third and forth argument note the direction of the dirac
        fnm = conj(AKsh(N, [], 0, 0)).';
    case 'dipole'
        fnm = zeros((N+1)^2,1);
        fnm(3) = 1;
end

clear nn Pm1 Pp1

% The function returns the correct coefficients, but one coefficient in [1]
% first edition is wrong! (See erratum to the book)
[gnm, D] = AKshRotate(fnm, rotation_angles);
% the rotation can also be applied using matrix multiplication (eq. (1.78), in [1])
% gnm = D * fnm;

% get the functions [1], eq. (1.40, 3.26)
F = Ynm * fnm;
G = Ynm * gnm;

% plot
AKf(20,10)
subtightplot(1,2,1, [0 0])
subtightplot(1,2,1, [0 0])
AKp(real(F)', 'x7', 'az', g(:,1), 'el', g(:,2), 'cb', 0, 'sph_proc', 'tri', 'coord', 2, 'hp_view', [180 0])
subtightplot(1,2,2, [0 0])
AKp(real(G)', 'x7', 'az', g(:,1), 'el', g(:,2), 'cb', 0, 'sph_proc', 'tri', 'coord', 2, 'hp_view', [180 0])


%% --------- 4. calculate the energy using spherical harmonics coefficients

% The energy of a spherical function f(phi, theta) can be calculated using
% the SH coefficients and Parseval's theorem (Eq. (1.43 in [1])

% In this case we use this to calcualte the diffuse field transfer function
% of an HRTF dataset and check the ennergy across SH orders
AKdependencies('FABIAN')
load('FABIAN_HRIR_modeled_HATO_0.mat', 'SH')

[e_nm, ~, e_N, avg_N] = AKshEnergy(SH.coeffLeft);

AKf(20,30)

% plot diffuse field transfer function
subplot(3,1,1)
    semilogx(SH.f, 10*log10(e_nm))
    axis([SH.f(2) 20e3 -15 15])
    title  'Diffuse field transfer function of FABIAN''s left ear'
    xlabel 'f in Hz'
    ylabel 'Magnitude spectrum in dB'
    grid on

% plot energy across SH orders for selected frequencies
subplot(3,1,2)
    f    = [1e3 2e3 4e3 8e3 16e3];
    f_id = round(f / (44100/256)) + 1;
    plot(0:SH.order, 10*log10(e_N(:,f_id)))
    hold on
    plot(0:SH.order, 10*log10(avg_N), 'k', 'linewidth', 2)
    legend('1 kHz', '2 kHz', '4 kHz', '8 kHz', '16 kHz', 'avg. across all freq.', 'location', 'SouthWest')
    title  'Energy across SH order for selected frequencies'
    xlabel 'SH order'
    ylabel 'energy in dB'
    grid on
    
% energy loss due to order truncation
subplot(3,1,3)
    trunc_error = cumsum( flipud( [e_N(:,f_id) avg_N] ) );
    trunc_error = AKm( flipud(trunc_error), max(trunc_error), '/');
    plot(0:SH.order, trunc_error(:,1:end-1))
    hold on
    plot(0:SH.order, trunc_error(:,end), 'k', 'linewidth', 2)
    title  'Energy decay across SH order for selected frequencies'
    xlabel 'SH order'
    ylabel 'energy in percent'
    grid on
