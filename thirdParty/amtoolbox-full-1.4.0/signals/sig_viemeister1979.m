function intervals = sig_viemeister79(trackvar,modfreq,bandwidth)
  
%   #Author : Peter L. Sondergaard  
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_viemeister1979.php

  
% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


fs=44100;
  
lowercutoff = 4000 - bandwidth;
uppercutoff = 4000;

siglen=22050;

nintervals=3;

intervals=randn(siglen,nintervals);

% Transform to frequency domain, and bandpass by cutting.

fint=fftreal(intervals);

cutvec=zeros(siglen/2+1,nintervals);
sp=round(lowercutoff/2);
up=round(uppercutoff/2);
cutvec(sp:up,:)=ones(up-sp+1,nintervals);

fint=fint.*cutvec;

intervals=ifftreal(fint,siglen);

% Amplitude modulate the target
intervals(:,1)=intervals(:,1).*(1 + (10^(trackvar/20)*sin(2*pi*modfreq/fs*(0:siglen-1)')));

% Apply starting and ending ramps.
intervals=rampsignal(intervals,round(siglen*.05));

% Adjust level for presentation
intervals=setdbspl(intervals,70);
  

