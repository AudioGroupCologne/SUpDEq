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
% S/W/G Sampled Wave Generator Wrapper R13-0306
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% [fftdata, kr] = sofia_swg(r, gridData, ac, FS, NFFT, ...
%                                       AZ, EL, Nlim, t, c, wavetype, ds)
% ------------------------------------------------------------------------
% fftData  Complex sound pressures                   [(N+1)^2 x NFFT]
% kr       kr-Vector
%          Can also be a matrix [krm; krs] for rigid sphere configurations:
%           [1,:] => krm referring to the Microphone Radius
%           [2,:] => krs referring to the Sphere/Microphone2 Radius
% ------------------------------------------------------------------------
% r        Microphone Radius
%          Can also be a vector for rigid sphere configurations:
%           [1,1] => rm  Microphone Radius
%           [2,1] => rs  Sphere Radius (Scatterer)
% gridData Quadrature grid                           [default LEB110]
%           Columns : Position Number 1...M
%           Rows    : [AZ EL Weight]
%           Angles AZ, EL in [RAD]
% ac       Array Configuration
%           0  Open Sphere with p Transducers
%           1  Open Sphere with pGrad Transducers
%           2  Rigid Sphere with p Transducers
%           3  Rigid Sphere with pGrad Transducers (Thx to Nils Peters!)
%           4  Dual Open Sphere with p Transducers (Thx to Nils Peters!)
% FS       Sampling Frequency
% NFFT     FFT Order (Number of bins) should be 2^x, x=1,2,3,...
% AZ       Azimuth   angle in [DEG] 0-2pi
% EL       Elevation angle in [DEG] 0-pi
% Nlim     Internal generator transform order limit
% t        Time Delay in s
% c        Speed of sound in [m/s] (Default: 343m/s)
% wavetype Type of the Wave. 0: Plane Wave (default) 1: Spherical Wave
% ds       Distance of the source in [m] (For wavetype = 1 only)
%          Warning: If NFFT is smaller than the time the wavefront 
%          needs to travel from the source to the array, the impulse 
%          response will by cyclically shifted (cyclic convolution).
% 
% This file is a wrapper generating the complex pressures at the
% positions given in 'gridData' for a full spectrum 0-FS/2 Hz (NFFT Bins)
% wave impinging to an array. The wrapper involves the W/G/C wave
% generator core and the I/T/C spatial transform core.
% 
% S/W/G emulates discrete sampling. You can observe alias artifacts.
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


function [fftData, kr] = sofia_swg(r, gridData, ac, FS, NFFT, AZ, EL, Nlim, t, c, wavetype, ds)

disp('SOFiA S/W/G Sampled Wave Generator Wrapper R13-0306');

if nargin < 12
    ds = 1;       %Default 1m
end

if nargin < 11
    wavetype = 0; %Default 0: Plane Wave
end

if nargin < 10
    c=343;        %Default 343m/s
end

if nargin < 9
    t=0;          %Default no timeshift
end

if nargin < 8
    Nlim = 120;   %Internal modal resolution limit of the wavefield generator
end

if nargin < 7
    EL = pi/2;    %Default 90°
end

if nargin < 6
    AZ = 0;       %Default 0°
end

if nargin < 5
    NFFT = 512;   %Default 512 Blocks
end

if nargin < 4
    FS = 48000;   %Default 48KHz
end

if nargin < 3
    ac = 0;       %Default Open Sphere
end

if nargin < 2
    gridData = sofia_lebedev(110, 0);  %EXAMPLE GRID LEB110, No Plot
end

if nargin == 0
    r = 0.1;
end

pause(1e-12);

if size(r, 2) == 2
   r(1,1) = r(1,1);
   r(2,1) = r(1,2);
   r(:,2) = [];
end

if size(r, 1) == 1
    kr = linspace(0,r*pi*FS/c,NFFT/2+1);
    krRef = kr;
else  
    kr(1,:) = linspace(0,r(1,1)*pi*FS/c,NFFT/2+1);
    kr(2,:) = linspace(0,r(2,1)*pi*FS/c,NFFT/2+1);
    
    if r(2,1) > r(1,1)
       krRef = kr(2,:);
    else
       krRef = kr(1,:);        
    end    
end

minOrderLim                     = 70; %Default minimum
rqOrders                        = ceil(krRef*2);
maxReqOrder                     = max(rqOrders);
rqOrders(rqOrders<=minOrderLim) = minOrderLim;
rqOrders(rqOrders>Nlim)         = Nlim;
Ng                              = max(rqOrders);
Pnm                             = zeros((Ng+1)^2,(NFFT/2+1));


index  = 1;
ctr    = -1;

if (maxReqOrder > Ng)
    disp(' ');
    disp(['   -> Warning: The requested wave needs a generator order of ', num2str(maxReqOrder),' to']);
    disp(['      deliver a fully valid result for the current configuration.']);
    disp(['      The generator is actually limited to a maximum order of ',num2str(Ng),'.']);
    disp(' ');
else
    if Ng == minOrderLim
        disp(['      Full spectrum generator order: ', num2str(Ng)]);
    else
        disp(['      Segmented generator orders: ', num2str(minOrderLim),' to ',num2str(Ng)]);
    end
end

while(1)  %Segmentation
    
    if mod(ctr,36) == 0 && ctr>0
        fprintf('\n      ');
    end
    fOrders = find(rqOrders == rqOrders(index));
    try
        Pnm = Pnm+sofia_wgc(Ng, r, ac, FS, NFFT, AZ, EL, t, c, wavetype, ds, min(fOrders), max(fOrders), rqOrders(index));
    catch ME
        if strcmp(ME.identifier,'MATLAB:unexpectedCPPexception')
            error('SOFiA:WGC:WaveUnstable','SOFiA W/G/C: Result is unstable. Change simulation parameters (r, FS, Nlim).');
        else
            error(ME.identifier,ME.message);
        end
        fftData = zeros(size(gridData,1),(NFFT/2+1));
        return
    end
    index = index+size(fOrders,2);
    if index > size(rqOrders,2)
        break;
    end
    ctr=ctr+1;
    pause(1e-12)
end

pause(1e-12);

try
    fftData = sofia_itc(Pnm, gridData);
catch ME
    if strcmp(ME.identifier,'MATLAB:unexpectedCPPexception')
        error('SOFiA:WGC:WaveUnstable','SOFiA W/G/C: Result is unstable. Change simulation parameters (r, FS, Nlim).');
    else
        error(ME.identifier,ME.message);
    end
    fftData = zeros(size(gridData,1),(NFFT/2+1));
    return
end

