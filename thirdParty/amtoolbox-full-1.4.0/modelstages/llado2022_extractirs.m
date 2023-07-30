function [IRs] = llado2022_extractirs(sofaFileName,lat_angles,fs)
%LLADO2022_EXTRACTIRS extract impulse response signals for the given lateral angles from a sofa file
%
%   Input parameters:
%     sofaFileName   : filename of the HRIR set in sofa format
%     lat_angles     : vector of lateral angles to extract
%     fs             : sample frequency
%
%   Output parameters:
%     IRs            : Impulse response according to the following matrix
%                      dimensions: direction x time x channel/ear
%
%   LLADO2022_EXTRACTIRS extracts the impulse responses as a SOFA file
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/llado2022_extractirs.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB M - Communication Systems
%   #Author: Pedro Llado (2022)
%   #Author: Petteri Hyvärinen (2022)
%   #Author: Ville Pulkki (2022)

% This file is licensed unter the GNU General Public License (GPL) either
% version 3 of the license, or any later version as published by the Free Software
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and
% at <https://www.gnu.org/licenses/gpl-3.0.html>.
% You can redistribute this file and/or modify it under the terms of the GPLv3.
% This file is distributed without any warranty; without even the implied warranty
% of merchantability or fitness for a particular purpose.

hrirname = char(sofaFileName);
SOFAobj = amt_load('llado2022',hrirname);
if (fs ~= SOFAobj.Data.SamplingRate)
    amt_disp('Sample Frequency mismatch');
end

for lat_id = 1:length(lat_angles)
    lat = lat_angles(lat_id);
    idx=find(round(SOFAobj.SourcePosition(:,1))==lat & SOFAobj.SourcePosition(:,2)==0,1);
    IRs(lat_id,:,:)=squeeze(SOFAobj.Data.IR(idx,:,:));
end
IRs = permute(IRs,[1,3,2]);
end



