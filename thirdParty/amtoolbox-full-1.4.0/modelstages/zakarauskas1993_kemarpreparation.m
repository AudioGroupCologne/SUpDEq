function hM = zakarauskas1993_kemarpreparation(varargin)
%ZAKARAUSKAS1993_KEMARPREPARATION Kemar HRTF preprocessing
%
%   Usage: hM = zakarauskas1993_kemarpreparation;
%
%   Preparation of KEMAR data to validate Zakarauskas et al. (1993)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/zakarauskas1993_kemarpreparation.php


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



amt_load('zakarauskas1993','hrtf_M_KEMAR normal pinna 22kHz');
hM=double(hM);
h=zeros(256,710,2);

for jj=1:size(hM,2)
    for ch=1:2
        h(:,jj,ch)=resample(hM(:,jj,ch),11,24);
    end
end
hM=h;
hM=hrtf2dtf(hM);  
hM=single(hM);


