function [BW, QERB]= erbest(impulse_ch, fsig, Fs)
%ERBEST estimate ERB from impulse response
%
%   ERBEST accepts an impulse response as an input. Based on its 
%   FFT and calculates its ERB (Equivalent Rectangular Bandwidth) 
%   and quality factor (QERB).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/erbest.php


%   #Author: Dick Lyon (2011)
%   #Author: Clara Hollomey (2021): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
 
le=length(impulse_ch);
F=(2*abs(fft(impulse_ch))/le).^2;
[M,Mn]=max(F(1:floor(le/2)));
E=sum(F(1:floor(le/2)));
BW=(E/M)*Fs/le;
QERB=fsig/BW;


