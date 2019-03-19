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
% Spherical Data Import R13-0306
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% timeData = sofia_readData(impulseResponses, FS, temperature,...
%                                              [downSample], [normalize])
% ------------------------------------------------------------------------     
% timeData            Struct with fields:
% 
%                     .impulseResponses     [N:Channels M:Samples]          
%                     .FS
%                     .downSample
%                     .averageAirTemp       Temperature in DEG 
%                     .irOverlay            Plot this for a good total 
%                                           overview of the dataset.
% ------------------------------------------------------------------------
% impulseResponses   Impulse response data [M:Channels N:Samples]
%                    !The audio data must be in in columns as it is  
%                    usual in MATLAB. But SOFiA works on audio data 
%                    in rows as it is usual in audio editors etc.
%                    Therefore the output IR matrix is flipped.
% 
%                    WARNING: Be sure the IRs come in the same order  
%                    as the sample grid points are organized!
% 
% FS                 Sampling frequency of the source material in 1/s
% 
% radius             Sensor sphere radius in meters
% 
% temperature        Average air temperature during the capture
%                    This is important for F/D/T to calculate the
%                    kr-vector. (This value is simply passed thru).
% 
% downSample         Downsampling factor  [default = 1: No downsampling]
%                    Downsampling is done using DECIMATE and a FIR low 
%                    pass filter of order 30. See MATLAB documentation 
%                    for more information. 
%                    !!! MATLAB Signal Processing Library required
% 
% normalize          Normalize flag 1:on, 0:off         [default = 1: on]       
%                    Normalizes the impulse responses with respect to the 
%                    absolute maximum value within the complete dataset.
% 
% directory          VSA dataset directory. If not defined a user dialog
%                    opens to pick a directory.
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


function timeData = sofia_readData(impulseResponses, FS, temperature, downSample, normalize)
    
disp('SOFiA Spherical Data Import R13-0306');

if nargin < 2
   error('Input arguments missing: timeData = sofia_mergeArrayData(impulseResponses, FS, [temperature], [downSample], [normalize])');
end

if nargin < 3
   temperature = 20; %Default 
end

if nargin < 4
   downSample = 1; %Default 
end

if nargin < 5
   normalize = 1; %Default 
end
  
if downSample > 1
   impulseResponses=decimate(cast(impulseResponses,'double'),downSample ,'FIR');
end

if normalize == 1
   impulseResponses = impulseResponses./max(max(abs(impulseResponses)));     
end

%Combine struct

timeData.impulseResponses = impulseResponses';
timeData.FS               = FS;
timeData.downSample       = downSample;
timeData.averageAirTemp   = temperature;
timeData.irOverlay        = sum(timeData.impulseResponses,1);
timeData.irOverlay        = abs(timeData.irOverlay/max(abs(timeData.irOverlay)));

disp([num2str(size(impulseResponses,2)),' impulse responses merged.']); 

if size(impulseResponses,2) > size(impulseResponses,1)
   disp(' ');
   warning('Impulse responses must be delivered as [M:Channels N:Samples]. The dimensions of your IRs are suspicious...')
end

fprintf('\n');
