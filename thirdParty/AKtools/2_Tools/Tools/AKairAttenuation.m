% [alpha, m, accuracy] = AKairAttenuation(f, h_r, T, p_a, doPlot)
% calculates the air attenuation in dB per meter according to [1]
%
% for example
% AKairAttenuation(1000)
% returns the air attenuation due to a travel distance of 1 m at a
% frequency of 1 kHz. To obtain the attanuation at 10 m multiply the result
% by 10.
% Empty bracktes can be passed to use a default parameter, e.g.
% AKairAttenuation(1000, [], 30)
% uses the default relative humidity, and a temperature of 30 degree.
%
%
% I N P U T
% f      - frequncy in Hz (scalar or vector)
% h_r    - relative humidity percent (default = 50)
% T      - temperature in degree Celsius (default = 20)
% p_a    - atmospheric pressure in Pa (default = 101325)
% doPlot - plot the valid range of atmospheric pressure vs. frequency
%          according to Section 7 in [1]
%
%
% O U T P U T
% alpha    - air attenuation per meter [dB] as specified by ISO 9613-1
% m        - energetic air attenuation per meter [linear] as needed for the
%            estimation of the reverberation time (term 4mV) as specified
%            by [2, Table 6.1, Eq. (5.24), and (5.25)]
% accuracy - string that gives the accuracy according to Section 7 in [1]
%
%
% [1] ISO 9613 "Attenuation of sound during propagation outdoors. Part 1:
%     Calculation of the absorption of sound by the atmosphere."
% [2] Heinrich Kuttruff, Room Acoustics, 5th ed., Spoon Press, Oxon, 2009.
%
% 2018/02 - fabian.brinkmann@tu-berlin.de

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

function [alpha, m, accuracy] = AKairAttenuation(f, h_r, T, p_a, doPlot)

% --- default parameter ----
if ~exist('p_a', 'var')
    p_a = 101325; % Pa
end
if isempty(p_a)
    p_a = 101325; % hPa
end
if ~exist('T', 'var')
    T = 20; % Celcius
end
if isempty(T)
    T = 20; % Celcius
end
if ~exist('h_r', 'var')
    h_r = 50; % percent
end
if isempty(h_r)
    h_r = 50; % percent
end
if ~exist('doPlot', 'var')
    doPlot = false;
end

% ---- reference values ----
% reference temperature in Kelvin
T_0  = 293.15;
% tripple point isotherm temperature in Kelvin
T_01 = 273.16;
% reference ambience atmospheric pressure in kPa
p_r  = 101325;

% ---- parameter conversions ----
% T in Kelvin
T = T + 273.15;

% relative humidity (h_r) to concentration of water vapour in percent (h)
C     = -6.8346 * (T_01/T)^1.261 + 4.6151;  % Eq. (B3)
p_sat = 10^C * p_r;                         % Eq. (B2)
h     = h_r *  p_sat / p_a;                 % Eq. (B1)

% ---- calculate the air attenuation ----
% Eq. (3)
f_rO = p_a/p_r * (24 + 4.04e4 * h * (.02 + h)/(.391 + h) );

% Eq. (4)
f_rN = p_a/p_r * (T/T_0)^(-1/2) * ...
        ( ...
            9 + 280*h * exp(-4.17 * ( (T/T_0)^(-1/3) - 1 ) ) ...
        );

% Eq. (5) in neper/m
alpha = f.^2 .* ...
    ( ...
        (1.84e-11 * (p_r/p_a) * (T/T_0)^(1/2)) + ...
        (T/T_0)^(-5/2)                         * ...
            ( ...
                .01275 * exp(-2239.1/T)     * ...
                ( f_rO + (f.^2/f_rO) ).^-1  + ...
                .1068 * exp(-3352/T)        * ...
                ( f_rN + (f.^2/f_rN) ).^-1    ...
            ) ...
    );

% convert to dB/m (factor 8.686)
alpha = alpha * 20*log10( exp(1) );

% convert to linear energetic value (thus the ten times logarithm)
m = alpha / ( 10*log10( exp(1) ) );

%  ---- estimate the accuracy according to section 7 in [1] ----
if nargout == 3
    
    % related humidity
    if     h >= .05 && h <= 5           && ...
           T >= 253.15 && T <= 323.15   && ...
           p_a < 200000
   
        accuracy = '+/- 10 percent';
        
    elseif h >= .005 && h <= 5           && ...
           T >= 253.15 && T <= 323.15    && ...
           p_a < 200000
       
        accuracy = '+/- 20 percent';
    else
        accuracy = '+/- 50 percent';
    end
    
end

% --- plot valid atmospheric pressure vs. frequency ----
% according to section 7 in [1]
if doPlot
    f_plot = logspace(1, log10(20e3), 100);
    p_1    = f_plot / 40;
    p_2    = f_plot*1e-6;
    
    AKf(20,10)
    semilogx(f_plot, [p_1; p_2], 'k', 'LineWidth', 2)
    legend('low', 'high')
    axis([10 20e3 -.5 5.5])
    set(gca, 'YTick', 0:5)
    AKfrequencyTicks
    title  'valid atmospheric pressure range vs. frequency'
    xlabel 'frequency / Hz'
    ylabel 'pressure / 1013.25 hPa'
    grid on
end