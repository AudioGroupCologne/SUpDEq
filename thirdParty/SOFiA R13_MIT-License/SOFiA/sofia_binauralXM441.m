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
% BinauralXM441 - Binaural Synthesis to Miro R13-0306
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
%
% miroObj = sofia_binauralXM441(Pnm, dn, miroObj, normalize, cGridType,... 
%                                                   hpfOn, hpTaps, ncGap)
% ------------------------------------------------------------------------     
% miroObj            Miro output object
% ------------------------------------------------------------------------              
% Pnm                Spatial Fourier Coefficients
% 
% dn                 Modal Array Filters from SOFiA M/F
%                   
% miroObj            Prototype Miro input object (e.g. array set), if []
%                    a new empty prototype is created.
%
% normalize          Normalize final BRIR set: false (#default) / true
%                  
% [cGridType]        Composite grid type 
%                    0 Lebedev quadrature (#default)
%                    1 Gauss-Legendre quadrature
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
% icients Pnm like binauralX(), but binauralXM() directly writes into a 
% Miro BRIR object instead of retuning raw binaural data. A prototype Miro
% object can be passed in order to maintain previous meta information.
% The returned object provides binaural BRIR data (1° azimuth resolution) 
% that can e.g. be used for dynamic binaural synthesis.
%
% For further information please refer to binauralX().
%
% Dependencies: - sofia_binauralX 
%               - sofia_gauss 
%               - sofia_tdt
%               - Miro class definition: "miro.m" 
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


function miroObj = sofia_binauralXM441(Pnm, dn, miroObj, normalize, cGridType, hpfOn, hpTaps, ncGap)

if nargin < 8
    ncGap = 16; 
end

if nargin < 7
    hpTaps = 2048; 
end

if nargin < 6
    hpfOn = true; 
end

if nargin < 5
    cGridType = 0;
end

if nargin < 4
    normalize = false;
end

if nargin < 3
    miroObj = miro();
end

if isempty(miroObj.name)
    miroObj.name = 'BinX';
else
    miroObj.name = [miroObj.name,'_BinX'];
end

if isempty(miroObj.date)
    miroObj.date = datestr(now);
end

if isempty(miroObj.miroVersion)
    miroObj.miroVersion = 1;
end

if ~isempty(miroObj.fs) && miroObj.fs ~= 44100
       warning('Sampling rate of the source Miro object is not 44100!');
end

miroObj.fs    = 44100;
miroObj.type  = 'BRIR';
miroObj.chOne = 'Left Ear';
miroObj.chTwo = 'Right Ear';
miroObj.nIr   = 360;
miroObj.positionReference = 'Head Rotation';

binauralGrid       = sofia_gauss(360,1,0);
miroObj.azimuth    = binauralGrid(:,1)';
miroObj.elevation  = binauralGrid(:,2)';
miroObj.quadWeight = binauralGrid(:,3)';
miroObj.quadGrid   = 'Gauss-Leg. 360SP (1E/360A)';
miroObj.scatterer  = true; 
miroObj.radius     = 0.0875; 

if isempty(miroObj.comments) || strcmp(miroObj.comments,'-')
    miroObj.comments = ['BinX: ', num2str(ncGap) ,' samples of ncGap latency added.'];
else
    miroObj.comments = [miroObj.comments,'_BinX: ', num2str(ncGap) ,' samples of ncGap latency added.'];
end

if isempty(miroObj.postProcessing) || strcmp(miroObj.postProcessing,'-')
    miroObj.postProcessing = ['SOFiA BinauralX - Binaural Synthesis'];
else
    miroObj.postProcessing = [miroObj.postProcessing,'_ SOFiA BinauralX - Binaural Synthesis'];
end

[BinauralL, BinauralR] = sofia_binauralX(Pnm, dn, miroObj.azimuth, cGridType, hpfOn, hpTaps, ncGap);

miroObj.irChOne = sofia_tdt(BinauralL)';
miroObj.irChTwo = sofia_tdt(BinauralR)';

if normalize
    norMax(1) = max(max(miroObj.irChOne));
    norMax(2) = max(max(miroObj.irChTwo));
    miroObj.normalization = 0.99/max(norMax);
    miroObj.irChOne = miroObj.irChOne * miroObj.normalization;
    miroObj.irChTwo = miroObj.irChTwo * miroObj.normalization;    
else 
    miroObj.normalization = 1;
end

miroObj.headWin = round(ncGap/2); 
miroObj.tailWin = round(1/8*size(miroObj.irChOne,1)); 

if ~isempty(miroObj.irCenter)
    miroObj.irCenter = [zeros(ncGap,1); miroObj.irCenter];
    if size(miroObj.irCenter,1) < size(miroObj.irChOne,1)
        zeroFiller = size(miroObj.irChOne,1) - size(miroObj.irCenter,1); 
        miroObj.irCenter = [miroObj.irCenter; zeros(zeroFiller,1)];
    elseif size(miroObj.irCenter,1) > size(miroObj.irChOne,1)
        miroObj.irCenter = miroObj.irCenter(1:size(miroObj.irChOne,1));
    end
end

miroObj.taps       = size(miroObj.irChOne,1);
miroObj.returnTaps = size(miroObj.irChOne,1);