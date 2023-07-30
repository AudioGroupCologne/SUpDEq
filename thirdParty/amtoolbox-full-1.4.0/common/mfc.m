function mfc=mfc(fc)
%MFC Generate all possible modulation frequencies for a given fc
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/mfc.php


%   #Author: Peter L. Soendergaard (2011)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
s=[0,5,10*(5/3).^(0:9)];

Q = 2;
bw = 5;
ex=(1+1/(2*Q))/(1-1/(2*Q));

startmf = 5;

umf = min(fc.*0.25, 1000);  

tmp = fix((min(umf,10) - startmf)/bw);
tmp = 0:tmp;
mfc = startmf + 5*tmp;
tmp2 = (mfc(end)+bw/2)/(1-1/(2*Q));
tmp = fix(log(umf/tmp2)/log(ex));
tmp = 0:tmp;
tmp = ex.^tmp;
mfc=[0 mfc tmp2*tmp];

%OLDFORMAT


