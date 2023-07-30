function [dtf,ctf]=hrtf2dtf(hrtf,gate,f1,f2,fs)
%   HRTF2DTF extracts dtf [and ctf] out of hrtf data
%
%   Usage:        [dtf,ctf]=hrtf2dtf(hrtf)
%                 [dtf,ctf]=hrtf2dtf(hrtf,gate)
%                 [dtf,ctf]=hrtf2dtf(hrtf,gate,f1,f2,fs)
%
%   Input parameters:
%     hrtf:     complete hrtf data in ARI format (hM)
%     gate:     optional string for defining bounds of ctf division
%               'full' -> no bounds
%               'bounds' -> define bounds in f1 and f2
%               otherwise 60dB bounds will be calculated
%     f1:       lower frequency bound 
%     f2:       upper frequency bound
%     fs:       sampling frequency;                 default: 48kHz
%
%   Output parameters:
%     dtf:      directional transfer function
%     ctf:      common transfer function (similar for all source positions)
%
%   Derives the directional transfer function from a set of HRTFs.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/hrtf2dtf.php


%   #Author: Robert Baumgartner (2010)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% default settings
if ~exist('gate','var')
    gate='auto';
end
if ~exist('fs','var')
    fs=48000;
end

% ctf calculation
n=size(hrtf,1);
hrtff=fft(hrtf);
ctff=mean(20*log10(hrtff),2);

switch gate
    case 'full'
        idx=1:round(size(ctff,1)/2);
    case 'bounds'
        id1=round(n*f1/fs);
        id2=round(n*f2/fs);
        idx=id1:id2;
    case 'auto'
        idx=find(mean(abs(ctff),3) <= min(mean(abs(ctff),3))+60);
        idx=idx(1:round(length(idx)/2)); % positiv frequency part
end

ctff=10.^(ctff(:,:,:)/20);

% minimal phase
for ch=1:size(ctff,3)
      % decompose signal
    amp=abs(squeeze(ctff(:,1,ch)));
    anu=-imag(hilbert(log(amp))); % minimal phase
    an=anu-round(anu/2/pi)*2*pi;  % wrap around +/-pi: wrap(x)=x-round(x/2/pi)*2*pi
    ctff(:,1,ch)=amp.*exp(1i*an);
end

% extracting positiv frequency part
ctff=ctff(1:round(n/2),:,:);
hrtff=hrtff(1:round(n/2),:,:);

% dtf calculation
dtff=hrtff;
for jj=1:length(idx)
    nn=idx(jj);
    for ii=1:size(dtff,2)
        for ch=1:2
            dtff(nn,ii,ch)=hrtff(nn,ii,ch)./ctff(nn,1,ch);
        end
    end
end

ltfatstart;
dtf=ifftreal(dtff,n);
ctf=ifftreal(ctff,n);
if size(ctf,1)>=240
    ctf=ctf(1:240,:,:); % rect windowing
end
end


