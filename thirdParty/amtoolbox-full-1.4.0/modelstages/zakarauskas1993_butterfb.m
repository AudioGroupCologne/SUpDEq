function [ y1n,y2n,fn ] = zakarauskas1993_butterfb( x1,x2,q10,fstart,fend,fs,space )
%ZAKARAUSKAS1993_BUTTERFB 4th-order butterworth filter bank
%
%   Usage:        [y1n,y2n,fn]= zakarauskas1993_butterfb( x1,x2 )
%                 [y1n,y2n,fn]= zakarauskas1993_butterfb( x1,x2,q10,fstart,fend,fs,space )
%
%   Input parameters:
%     x1      : impulse responses
%     x2      : impulse responses
%     q10     : relativ bandwidth of filter bands;  default... 10
%     fstart  : start frequency; minimum: 0,5kHz;   default... 2kHz
%     fend    : end frequency;                      default... 16kHz
%     fs      : sampling frequency;                 default... 48kHz
%     space   : spacing factor of filter bands;     default... 1.03
%
%   Output parameters:
%     y1n:      dB-values of each filter band and input signal x1
%     y2n:      dB-values of each filter band and input signal x2
%     fn:       center frequencies of filter bands
%
%   Description:
%     4th-order butterworth filter bank in order to get approximately 
%     a bank of filters having parabolic shape according to Zakarauskas
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/zakarauskas1993_butterfb.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Robert Baumgartner (2010), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% default settings
if ~exist('q10','var')
    q10=10;
end
if ~exist('fstart','var')
    fstart=1000;
end
if ~exist('fs','var')
    fs=48000;
end
if ~exist('fend','var')
    fend=fs/2;
end
if ~exist('space','var')
    space=1.03;
end

if fs<=22050 && space <= 1.05
    hr=2;       %in order to keep F3dB2 under Nyquist frequency
elseif fs<=22050 && space <= 1.10 && space > 1.05
    hr=1;
else
    hr=0;
end
n=0:log(fend/fstart)/log(space)-hr;
fn=fstart*space.^n;
bwe=0.12/q10;
y1n=zeros(length(fn),size(x1,2),size(x1,3));
y2n=y1n;
for ii=1:length(fn)             % filter bands
    d=fdesign.bandpass('N,F3dB1,F3dB2',4,10^(log10(fn(ii))-bwe)*2/fs,10^(log10(fn(ii))+bwe)*2/fs);
    hd=design(d,'butter'); 
    for ch=1:2                  % channels
        for pos=1:size(y1n,2)   % positions
            y1temp=filter(hd,x1(:,pos,ch));
            y1n(ii,pos,ch)=10*log10(sum(y1temp.^2));
            y2temp=filter(hd,x2(:,pos,ch));
            y2n(ii,pos,ch)=10*log10(sum(y2temp.^2));
        end
    end
% uncommend the following 3 lines to plot the spectra of all filter bands
%             imp=zeros(100000,1);imp(1)=1;
%             y=filter(hd,imp);
%             spect(fft(y),fs)
end
end


