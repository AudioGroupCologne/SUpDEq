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
% T/D/T - Time Domain Transform R13-0306
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% timeDomainSignal = sofia_tdt(Y, [win], [resampleFactor], [minPhase])
% ------------------------------------------------------------------------     
% timeDomainSignal   Reconstructed Time Domain Signal
%                    Columns : Index / Channel: IR1, IR2, ... IRn
%                    Rows    : Impulse responses (time domain)
% ------------------------------------------------------------------------              
% Y                  Frequency domain FFT data for multiple channels
%                    Columns : Index / Channel
%                    Rows    : FFT data (frequency domain)
% 
% [win]              Window Signal tail [0...1] with a HANN window
%                    0    off (#default)
%                    0-1  window coverage (1 full, 0 off)
%                   
% [resampleFactor]   Optional resampling: Resampling factor 
%                    e.g. FS_target/FS_source
%                    Resampling is done using MATLAB RESAMPLE 
%                    (See MATLAB documentation for more details)
%                    ! Signal Processing Toolbox required  
%
% [minPhase]         Optional minimum phase reduction 
%                    0 off (#default)
%                    1 on
%                    ! Signal Processing Toolbox required 
%
% 
% This function recombines time domain signals for multiple channels from
% frequency domain data. It is made to work with half-sided spectrum FFT 
% data.  The impulse responses can be windowed.  The IFFT blocklength is 
% determined by the Y data itself:
% 
% Y should have a size [NumberOfChannels x ((2^n)/2)+1] with n=[1,2,3,...] 
% and the function returns [NumberOfChannels x resampleFactor*2^n] samples.
% 
% Dependencies: MATLAB Signal Processing Toolbox required for 
%               windowing and resampling.    
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

function timeDomainSignal = sofia_tdt(Y, win, resampleFactor, minPhase)

disp('SOFiA T/D/T - Time Domain Transform R13-0306');

if nargin == 0
   error('Arguments missing: timeDomainSignal = sofia_tdt(Y, [win], [resampleFactor])');
end

if nargin < 2
   win = 0;
end

if win > 1
   error('Argument win must be in the range 0 to 1.')
end

if nargin < 4
   minPhase = 0;
end

signalProcessingExists = 0;
toolboxes = ver;
for i = 1: size(toolboxes,2)
    if ~isempty(strfind(toolboxes(i).Name,'Signal Processing Toolbox'))
       signalProcessingExists = 1;   
       break
    end
end

if ~signalProcessingExists 
   disp('WARNING: Signal processing toolbox not found. Resampling and minimum phase disabled.') 
end

Y(:,1) = Y(:,2); 
y = [Y(:,:), conj(fliplr(Y(:,2:end-1)))];
y = real(ifft(y,[],2)); 

if signalProcessingExists && minPhase ~= 0    
    y    = [y, zeros(size(y))]';
    Y    = fft(y);
    Y(Y == 0) = 1e-21;
    img  = imag(hilbert(log(abs(Y))));
    y    = real(ifft(abs(Y) .* exp(-1i*img)));
    y    = y(1:end/2,:)';    
end

if win ~= 0   
    irLength = size(y,2);
    j = floor(irLength*win):irLength;
    winfkt = 0.5+0.5*cos(2*pi*(j-((irLength-1)/2))/(irLength-1)); 
    y(:,end-size(winfkt,2)+1:end) = y(:,end-size(winfkt,2)+1:end).*repmat(winfkt,size(y,1),1);
end

if nargin == 3 && signalProcessingExists 
    timeDomainSignal = resample(y',resampleFactor,1)';
else
    timeDomainSignal = y;   
end

