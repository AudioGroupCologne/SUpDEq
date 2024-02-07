function [dSRTs_out, newSRT] = joergensen2011_pctodsrt(PC_input,SNRs,conditions)
%JOERGENSEN2011_PCTODSRT calculates the SRT and change in SRT
% 
%   Usage:
%     [dSRTs_out, newSRT] = joergensen2011_pctodsrt(PC_input,SNRs,conditions)
%
%   Input parameters:
%     PC_input   :  Matrix with the mean percent correct for each processing
%                   condition and SNR. The first column of PC_input should always be the
%                   reference with no processing
%     conditions :  Vector with the processing conditions
%     SNRs       :  Vector with the SNRs used.
%
%   JOERGENSEN2011_PCTODSRT calculates the SRT and change in SRT from the simulated percent correct 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/joergensen2011_pctodsrt.php


%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Untrusted
%   #Requirements: M-Signal M-Stats
%   #Author: Søren Jørgensen (2011)
%   #Author: Peter L. Sondergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

 % ----------------  calculating dSRTs on mean psychofuncs.
for k = conditions
    SNRsLong = 0;
    yhat = 0;
    % ----------- Connecting the points with streight lines-----------
    for p = 1:length(SNRs)-1
        tmp =  polyfit(SNRs(p:p+1), PC_input(p:p+1,k)',1);
        tmpSNR =  SNRs(p):0.005:SNRs(p+1);
        yhat_tmp = polyval(tmp,tmpSNR);
        
        yhat = [yhat yhat_tmp];
        SNRsLong = [SNRsLong tmpSNR];
        
    end
    
    yhat2(:,k) = yhat(2:end);
    
    SNRsLong = SNRsLong(2:end);
    SRT = 50;
    hLine = linspace(SRT,SRT,length(yhat2(:,k)))';
    C = abs(yhat2(:,k)-hLine); %  Find minimum distance between horizontal line and estimated p_correct
    SNRIndex = find(C==min(C));
    if length(SNRIndex)>1
        SNRIndex = SNRIndex(1);
    end
    newSRT(k) =  SNRsLong(SNRIndex);
    
%     ----------- calculating change in SRT ---------------------
    dSRTs = abs(newSRT(1)-newSRT(k));
    if newSRT(k) > newSRT(1)
        dSRTs = dSRTs;
    else
        dSRTs = dSRTs*-1;
    end
    if min(C) > 1
        dSRTs = -4;
    end
    if  isempty(dSRTs)
        dSRTs = -4;
    end
    dSRTs_out(k) = dSRTs;
end


