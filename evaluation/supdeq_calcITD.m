%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [itd_ms, itd_ys] = supdeq_calcITD( hrir, fs, th, us, fc )
%
% This function calculates the ITD (Interaural Time Delay) of a HRIR based
% on onset detection 
%
% Output:
% itd_ms            - ITD of the HRIR(s) in ms 
% itd_ys            - ITD of the HRIR(s) in micro seconds 
%
% Input:
% hrir              - [N x 2] HRIR with N samples and 2 channels    
% fs                - Sampling Rate
%                     Default: 48000
% th                - Threshold (in dB) of the onset detection
%                     Default: -20
% us                - Upsampling factor for onset detection
%                     Default: 10
% fc                - Frequency in Hz at which an 8th order Butterworth 
%                     lowpass is applied to the hrir before onset detection
%                     Default: 3000
%
% Dependencies: AKTools
%
% Reference:
% A. Andreopoulou and B. F. G. Katz, ?Identification of perceptually relevant
% methods of inter-aural time difference estimation,? 
% J. Acoust. Soc. Am., vol. 142, no. 2, pp. 588?598, 2017.
%   
% (C) 2020 by JMA, Johannes M. Arend  
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [itd_ms, itd_ys] = supdeq_calcITD( hrir, fs, th, us, fc )

if nargin < 2 || isempty(fs)
    fs = 48000;
end

if nargin < 3 || isempty(th)
    th = -20;
end

if nargin < 4 || isempty(us)
    us = 10;
end

if nargin < 5 || isempty(fc)
    fc = 3000;
end

%Check type of array
if size(hrir,2) > size(hrir,1)
    hrir = hrir';
end

%%

onsetL = AKonsetDetect(hrir(:,1), us, th, 'rel', [fc fs]);
onsetR = AKonsetDetect(hrir(:,2), us, th, 'rel', [fc fs]);

itdSamples = onsetL-onsetR;
itd_ms = itdSamples/fs*1000;
itd_ys = itd_ms*1000;
    
end

