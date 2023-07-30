function sig = sig_schroeder1970(f0,n,C,fs,D)
%SIG_SCHROEDER1970 generates a Schroeder-phase harmonic complex tone
%   Usage: sig = sig_schroeder1970(f0,n,C,fs,D)
%
%   Input parameters:
%     f0 : fundamental frequency (Hz)
%     n : index vector specifying the contained harmonics
%     C : phase curvature, [-1,1]
%     fs : sampling rate (Hz)
%     D : signal duration (sec)
%
%   SIG_SCHROEDER1970 generates the Schroeder-signal with modified phase of the overtones
%
%   References:
%     M. R. Schroeder. Synthesis of low peak-factor signals and binary
%     sequences with low autocorrelation. IEEE Trans. Inf. Theory, 16:85--89,
%     1970.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_schroeder1970.php


%   #Author: Robert Baumgartner

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

t = 0:(1/fs):D;
phi = C*pi*n.*(n+1)/length(n);
f = f0*n;
sig = mean(repmat(f(end)./f(:),1,length(t)).*sin(2*pi*n'*f0*t + repmat(phi',1,length(t))),1);

end


