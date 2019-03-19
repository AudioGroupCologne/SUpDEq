% c = AKspeedOfSound(T, h_r, p_a, x_c)
%
% caluclates the speed of sound c [m/s] for a given temperature
% T in Celcius. If only the temperature is passed, this is done according
% to Eq. 2.17 in [1]: c = sqrt(kappa*R/M_mol*T_in_kalvin). Otherwise the
% approximation from Eq. (15) [2] is used. The approximation is valid for
% Temperature between 0 and 30 degree Celcius.
%
% I N P U T
% T   - temperature in degree Celcius
% h_r - relative humidity in percent
% p_a - atmospheric pressure in Pa (default = 101325)
% x_c - CO_2 concentration in air in percent  (defualt = 0.314)
%
%
% O U T P U T
% c   - speed of sound in m/s
%
%
% C O N S T A N T S
% M_mol       = 28.8e-3;          % molecular mass [g]
% kappa       = 1.4;              % adiabatic index
% R           = 8.314;            % gas constant [Nm/K]
%
%
% [1] Michael Möser (2007): "Technische Akustik", 7th edition, Eq. 2.17
% [2] O. Cramer (1993): "The variation of the specific heat ratio and the
%     speed of sound in air with temperature, pressure, humidity, and CO2
%     concentration'', J. Acoust. Soc. Am., 93(5): 2510-2516. 
%
%
% 03/2014  -  fabian.brinkmann@tu-berlin.de (initial dev.)
% 03/2018  -  fabian.brinkmann@tu-berlin.de (added calculation in
%             dependency of h_r, p_a, and x_c based on code from Stefan
%             Weinzierl)

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
function c = AKspeedOfSound(T, h_r, p_a, x_c)

M_mol  = 28.8e-3;          % molecular mass [g]
kappa  = 1.4;              % adiabatic index
R      = 8.314;            % gas constant [Nm/K]

if nargin == 1
    
    T = T+273.15;         % temperature in Kelvin
    
    % see Moeser, Michael (2007): Technische Akustik, 7th edition, Eq. 2.17
    c = sqrt(kappa*R/M_mol*T);
    
else
    
    % default parameter
    if ~exist('x_c', 'var')
        x_c = .314;
    end
    if ~exist('p_a', 'var')
        p_a = 101325; % Pa
    end
    
    x_c = x_c / 100;
    
    % Correction for wet air (valid for -50 < T in Celsius < 90)
    fw = 1.00519;
    
    % saturation vapour pressure of water with correction for wet air
    % according to: D. Sonntag, and D. Heinze (1982): Sättigungsdampfdruck-
    % und Sättigungsdampfdichtetafeln für Wasser und Eis. (1. Aufl.),
    % VEB Deutscher Verlag für Grundstoffindustrie
    % (Magnus-Formel)
    pws = fw * 611.213 * 10.^( (7.602*T) ./ (241.2+T) );
    
    % mixing ratio (mole fraction) of water vapor
    xw = 0.01*h_r * pws / p_a;
    
    % Coefficients according to [2]
    a0 = 331.5024;
    a1 = 0.603055;
    a2 = -0.000528;
    a3 = 51.471935;
    a4 = 0.1495874;
    a5 = -0.000782;
    a6 = -1.82e-7;
    a7 = 3.73e-8;
    a8 = -2.93e-10;
    a9 = -85.20931;
    a10 = -0.228525;
    a11 = 5.91e-5;
    a12 = -2.835149;
    a13 = -2.15e-13;
    a14 = 29.179762;
    a15 = 0.000486;
    
    % approximation for c according to [2]
    c1 = a0+a1*T+a2*T.^2;
    c2 = (a3+a4*T+a5*T.^2).*xw;
    c3 = (a6+a7*T+a8*T.^2).*p_a;
    c4 = (a9+a10*T+a11*T.^2).*x_c;
    c5 = a12*xw.^2+a13*p_a.^2+a14*x_c.^2+a15*xw.*p_a.*x_c;
    
    c = c1 + c2 + c3 + c4 + c5;  
    
    % check validity
    if T<0 || T > 30
        warning('AKspeedOfSound:Input', 'The approximation is only valid for 0<=T<=30')
    end
    if p_a<75e3 || T > 102e3
        warning('AKspeedOfSound:Input', 'The approximation is only valid for 75e3<=p_a<=102e3')
    end
    if x_c < 0 || x_c > .01
        warning('AKspeedOfSound:Input', 'The approximation is only valid for 0<=x_c<=1')
    end
    if xw > .06
        warning('AKspeedOfSound:Input', 'The approximation is only valid x_w<6 percent')
    end
    
end