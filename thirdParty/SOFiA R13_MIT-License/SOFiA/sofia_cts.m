% /// ASAR/MARA Research Group
%  
% Technology Arts Sciences TH Köln
% Technical University of Berlin
% Deutsche Telekom Laboratories
% University of Rostock
% WDR Westdeutscher Rundfunk
% IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis
% 
% C/T/S Center Transducer Simulation R13-0306
%
% Produces the corresponding center transducer signals for the wave
% generators. The center transducer signal is e.g. required for the
% SOFiA B/S/A BEMA Spatial Anti Aliasing module.
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% center = sofia_cts(FS, NFFT, AZ, EL, t, c)
% ------------------------------------------------------------------------
% center  Complex sound pressures (Center Transducer)  [1 x NFFT]
% ------------------------------------------------------------------------
% FS      Sampling rate [1/s]                       [default = 48000]
% NFFT    Number of FFT bins                        [default = 512]
% AZ      Azimuth angle in [RAD] 0-2pi              [default = 0]
% EL      Elevation angle in [RAD] 0-pi             [default = pi/2]
% t       Timeshift, for (FS/t) Samples             [default = 0]
% c       Speed of Sound [m/s]                      [default = 343]
%

% CONTACT AND LICENSE INFORMATION:
% 
% /// ASAR/MARA Research Group 
%  
%     [1] Technology Arts Sciences TH Köln
%     [2] Technical University of Berlin 
%     [3] Deutsche Telekom Laboratories 
%     [4] University of Rostock
%     [5] WDR Westdeutscher Rundfunk 
%     [6] IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis toolbox
% 
% Copyright 2011-2017 Benjamin Bernschütz et al.(§)  
% 
% Contact ------------------------------------
% Technology Arts Sciences TH Köln 
% Institute of Communications Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
% 
% phone       +49 221 8275 -2496 
% cell phone  +49 171 4176069 
% mail        rockzentrale 'at' me.com 
% --------------------------------------------
% 
% This file is part of the SOFiA sound field analysis toolbox
%
% Licence Type: MIT License
%
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
%
% (§) Christoph Pörschmann [1]     christoph.poerschmann 'at' th-koeln.de
%     Sascha Spors         [2,3,4] sascha.spors 'at' uni-rostock.de  
%     Stefan Weinzierl     [2]     stefan.weinzierl 'at' tu-berlin.de


function [center] = sofia_cts(FS, NFFT, AZ, EL, t, c)

if nargin < 6
    c=343;
end

if nargin < 5
    t=0;
end

if nargin < 4
    EL = pi/2;
end

if nargin < 3
    AZ = 0;
end

if nargin < 2
    NFFT = 512;
end

if nargin < 1
    FS = 48000;
end


if (FS*t >= NFFT)
    warning('Timeshift t > NFFT ! (Cyclic convolution)');
    disp(' ');
elseif (FS*t >= NFFT/2)
    warning('Timeshift t > NFFT/2 - Remember to guard headroom for filter responses.');
    disp(' ');
end

disp('SOFiA C/T/S - Center Transducer Simulation R13-0306');

EL = pi-EL;

f  = linspace(0,FS/2,NFFT/2+1);
k  = 2*pi*f/c;
w  = 2*pi*f;

fftBins = size(k,2);

center = zeros(1,fftBins);

[Xo,Yo,Zo]=sph2cart(0,0,0);
for bin=1:fftBins
    center(bin)=exp(1i*(k(bin)*Xo*sin(EL)*cos(AZ)+k(bin)...
        *Yo*sin(AZ)*sin(EL)+k(bin)*Zo*cos(EL)))*exp(-1i*w(bin)*t);
end


