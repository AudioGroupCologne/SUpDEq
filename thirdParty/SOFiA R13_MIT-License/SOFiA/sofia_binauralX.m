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
% BinauralX - Binaural Synthesis R13-0306
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
%
% [BinauralL, BinauralR] = sofia_binauralX(Pnm, dn, headRot, [cGridType], ...
%                          [hpfOn], [hpTaps], [ncGap])
% ------------------------------------------------------------------------     
% BinauralL/R        Binaural Frequency domain FFT data for multiple ch.
%                    Columns : Index / Channel
%                    Rows    : FFT data (frequency domain)
% ------------------------------------------------------------------------              
% Pnm                Spatial Fourier Coefficients
% 
% dn                 Modal Array Filters from SOFiA M/F
%                   
% headRot            Head rotation(s) in RAD
%                  
% [cGridType]        Composite grid type 
%                    0 Gauss-Legendre quadrature (#default)
%                    1 Lebedev quadrature 
%                    
% [hpfOn]            30Hz highpass filter 
%                    true  : HPF on (#default)
%                    false : HPF off
%
% [hpTaps]           Taps for the highpass filter, 2048 (#default)
%
% [ncGap]            Non-causality FIR gap (radial filters), 16 (#default)
%                     -> Entails latency of ncGap samples
%
%
% The function synhtesizes binaural signals from spatial fourier coeff- 
% icients Pnm. The binaural cues are generated using the HRTFs of a Neumann 
% KU100 artificial head. The function combines plane wave decomposition on 
% a matched order composite grid and spatial downsampling of high order 
% (N=35) interpolated HRTFs. The binaural signals can be directly calculated 
% for different azimutal head rotations in a single % run (e.g. to generate 
% full dynamic binaural synthesis datasets).
%
% ! IMPORTANT ADVICE:
%
%  I.  Only signals with a sampling rate of 48000Hz supported (HRTFs)
%  II. Initial NFFT must be at least twice the length of the IR
%    
%
% Dependencies: - sofia_pdc               
%               - sofia_itc
%               - sofia_lebedev or sofia_gauss
%               - SOFiA KU100 HRIR Set stored in: "sofia_HRIR.mat"
%               - MATLAB Signal Processing Toolbox for the HPF   
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


function [BinauralL, BinauralR] = sofia_binauralX(Pnm, dn, headRot, ...
                                  cGridType, hpfOn, hpTaps, ncGap)

if nargin < 7
    ncGap = 16; 
end

if nargin < 6
    hpTaps = 2048; 
end

if nargin < 5
    hpfOn = true; 
end

if nargin < 4
    cGridType = 0;
end

clc
fprintf('SOFiA BinauralX - Binaural Synthesis R13-0306\n\n');
fprintf(' -> Only signals with a sampling rate of 48000Hz supported!\n')
fprintf(' -> Initial NFFT must be at least twice the length of the IR!\n\n')
fprintf('Be patient - this may take a while...\n\n')
pause(1.5)

if ~exist('sofia_HRIR.mat','file')
    error('Binaural HRIR dataset not accessible! File: sofia_HRIR.mat');
else
    load sofia_HRIR
end

[compositeGrid, Npwd] = getCompositeGrid(Pnm, cGridType);

Y       = sofia_pdc(Npwd, compositeGrid(:,1:2), Pnm, dn);
weights = repmat(compositeGrid(:,3),1,size(Y,2));
Y       =  Y .* weights * (Npwd+1)^2;

[BinauralL, BinauralR] = binauralize(Hl_nm, Hr_nm, Y, compositeGrid, headRot, hpfOn, hpTaps, ncGap);



function [compositeGrid, Npwd] = getCompositeGrid(Pnm, cGridType)

Npwd = sqrt(size(Pnm,1))-1;

if cGridType
    
    lebNodes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
        350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, ...
        3074, 3470, 3890, 4334, 4802, 5294, 5810];
    
    lebOrders  = floor(sqrt(lebNodes/1.3)-1);
    matchOrder = abs(lebOrders-Npwd);
    orderIndex = find(matchOrder==min(matchOrder));
    
    if ~matchOrder(orderIndex)
        fprintf(['>>> composite grid: Lebedev, ', num2str(lebNodes(orderIndex)),' Nodes\n\n']);
        compositeGrid = sofia_lebedev(lebNodes(orderIndex),0);
        return
    else
        warning('No matching Lebedev composit grid available: Switched to Gauss grid.\n');
    end
end

fprintf(['>>> composite grid: Gauss, ', num2str(2*(Npwd+1)^2),' Nodes\n\n']);
compositeGrid = sofia_gauss(2*(Npwd+1), (Npwd+1),0);



function [BinauralL, BinauralR] = binauralize(Hl_nm, Hr_nm, Y, compositeGrid, headRot, hpfOn, hpTaps, ncGap)

NFFT = (2*(size(Y,2)-1));

Y = conj([Y(:,:), conj(fliplr(Y(:,2:end-1)))])';
y = real(ifft(Y));
y = circshift(y,ncGap); % cyclic move for ncGap

u      = 0:2*ncGap-1;   % window full block (radial filters)
winFkt = 0.5+0.5*cos(2*pi*(u-((2*ncGap-1)/2))/(2*ncGap-1));
winFkt = [winFkt(1:end/2), ones(1,NFFT/2), winFkt(end/2+1:end), zeros(1,NFFT/2-2*ncGap)]';
winFkt = repmat(winFkt,1,size(y,2));
y      = y.*winFkt;

if hpfOn %Recalc required NFFT
    NFFT = NFFT/2 + 2*ncGap + (size(Hl_nm,2)-1) + hpTaps; 
else
    NFFT = NFFT/2 + 2*ncGap + (size(Hl_nm,2)-1); %reserve 128taps only for HRTF 
end

if size(y,1) > NFFT
    y = y(1:NFFT,:);
end

Y = fft(y,NFFT);

BinauralL = zeros(length(headRot),NFFT);
BinauralR = zeros(length(headRot),NFFT);

fprintf('\n');

for hrPointer = 1:length(headRot)
    
    fprintf(['Head orientation ', num2str(mod(headRot(hrPointer)*180/pi,360)),'°\n']);
    hrtfRotAZ           = mod(compositeGrid(:,1)-headRot(hrPointer),2*pi);
    Hl                  = sofia_itc(Hl_nm, [hrtfRotAZ compositeGrid(:,2)]);
    Hr                  = sofia_itc(Hr_nm, [hrtfRotAZ compositeGrid(:,2)]);
    
    Hl   = conj([Hl(:,:), conj(fliplr(Hl(:,2:end-1)))])';
    Hr   = conj([Hr(:,:), conj(fliplr(Hr(:,2:end-1)))])';
    hl   = real(ifft(Hl));
    hr   = real(ifft(Hr));
    Hl   = fft(hl,NFFT);
    Hr   = fft(hr,NFFT);

    BinL = Y.*Hl;
    BinL = sum(BinL,2);
    BinauralL(hrPointer, :) = BinL;
    BinR = Y.*Hr;
    BinR = sum(BinR,2);
    BinauralR(hrPointer, :) = BinR;    
    
end

if hpfOn
    
    hpf       = 30;    
    fs        = 48000;
    hpf       = fir1(hpTaps,hpf/(fs/2),'high');
    HPF       = fft(hpf, 2*hpTaps);
    img       = imag(hilbert(log(abs(HPF))));
    hpf       = real(ifft(abs(HPF) .* exp(-1i*img)));
    hpf       = hpf(1:end/2);
    winLength = round(hpTaps/2);
    c         = 0:(round(winLength))-1;
    winFkt    = 0.5+0.5*cos(2*pi*(c-((winLength-1)/2))/(winLength-1));
    winFkt    = winFkt(end/2+1:end);
    winFkt    = [ones(1,hpTaps-size(winFkt,2)), winFkt];
    hpf       = hpf.*winFkt;
    HPF       = fft(hpf, NFFT);
    HPF       = repmat(HPF,size(BinauralL,1),1);
    BinauralL = BinauralL.*HPF;
    BinauralR = BinauralR.*HPF;
    
end

BinauralL = BinauralL(:,1:end/2+1);
BinauralR = BinauralR(:,1:end/2+1);


