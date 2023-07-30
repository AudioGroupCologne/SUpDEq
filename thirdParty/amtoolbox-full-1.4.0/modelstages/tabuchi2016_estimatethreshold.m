function [Kave, ThresEstVec] = tabuchi2016_estimatethreshold(kfunc, kvals)
%TABUCHI2016_ESTIMATETHRESHOLD Estimates the masked threshold              
%   Usage: [Kave, ThresEstVec] = tabuchi2016_estimatethreshold(kfunc, kvals)
%
%   Input parameters:
%     kfunc     : k function (k decision variable)
%     kvals     : k values (k decision variable)
%
%   Output parameters:
%     Kave         : k values averages
%     ThreshEstVec : estimated threshold in dB
%
%   estimates the threshold per masking condition
%
%   See also: tabuchi2016 exp_tabuchi2016
%
%   References:
%     H. Tabuchi, B. Laback, T. Necciari, and P. Majdak. The role of
%     compression in the simultaneous masker phase effect. The Journal of the
%     Acoustical Society of America, 140(4), 2016.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/tabuchi2016_estimatethreshold.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Verified
%   #Requirements: 
%   #Author: Hisaaki Tabuchi
%   #Author: Clara Hollomey (adaptations for AMT)


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

Cvec = -1.5:0.25:1;
GmaxVec = 0:70;
CondStr_vec = {'OffFreqPre', 'OnFreqPre'};
mlvl_vec = [60 90];

% load the k values
InK = kvals;
Kave = mean(InK(:));% the final Kave was 1.0261, not 1.06

% kfunc: dim1 for TargetLevel, dim2 for Cvec, dim3 for Gmax, dim4 for PrecCond, dim5 for MaskerLevel
InMat = kfunc;
% set up the matrix
ThresEstVec = NaN(size(InMat,2), size(InMat,3), size(InMat,4), size(InMat,5));

for mlvlLoop = 1:length(mlvl_vec)    
    for PrecCondLoop = 1:length(CondStr_vec)        
        for GmaxVecLoop = 1:length(GmaxVec)
            for CvecLoop = 1:length(Cvec)

                KvecMat_Interest = InMat(:, CvecLoop, GmaxVecLoop, PrecCondLoop, mlvlLoop);
                [minval,ind] = min(abs(KvecMat_Interest - Kave));  % find a taget dB for the minima of K(i.e., OutM+T - OutM) - Kave. 
                if minval > 0.1
                    error('minval is not appropriate.')
                end

                EstThresh_dB = (ind - 1)*0.1;   % in case the step size of target dB SPL (x_vec) is 0.1
                ThresEstVec(CvecLoop, GmaxVecLoop, PrecCondLoop, mlvlLoop) = EstThresh_dB; 

            end           
        end
    end    
end

end

