function [dSRTs_out, newSRT] = joergensen2011_PCtodSRT(PC_input,SNRs,conditions)
%joergensen2011_PCtodSRT: calculates the SRT and change in SRT from the simulated percent correct 
% 
% Usage: [dSRTs_out newSRT] = joergensen2011_PCtodSRT(PC_input,SNRs,conditions)
% PC_input   :  Matrix with the mean percent correct for each processing
%               condition and SNR. The first column of PC_input should always be the
%               reference with no processing
% conditions :  Vector with the processing conditions
% SNRs       :  Vector with the SNRs used.
% 
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/joergensen2011_PCtodSRT.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Was dSRTs_from_pc_mean by Soeren Joergensen august 2010
% roughly adapted to AMT by Piotr Majdak

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

