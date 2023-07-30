function sigma = baumgartner2016_comparisonprocess(tar,tem)
%BAUMGARTNER2016_COMPARISONPROCESS Comparison with direction-specific templates
%   Usage:     sigma = baumgartner2016_comparisonprocess(tar,tem)
%
%   Input parameters:
%     tar     : discharge rate profiles of target sounds (fields: tar.m for
%               magnitude and tar.sd for standard deviation)
%     tem     : discharge rate profiles of templates (fields as for tar)
%
%   Output parameters:
%     sigma   : internal distance metric
%
%   BAUMGARTNER2016_COMPARISONPROCESS(...) compares discharge rate profiles
%   on the basis of the quotient between rate difference and auditory nerve
%   variance (May and Huang, 1997; Reiss et al., 2011)
%
%   References:
%     B. J. May and A. Y. Huang. Spectral cues for sound localization in
%     cats: A model for discharge rate representation in the auditory nerve.
%     J. Acoust. Soc. Am., 101:2705--2719, 1997.
%     
%     L. A. J. Reiss, R. Ramachandran, and B. J. May. Effects of signal level
%     and background noise on spectral representations in the auditory nerve
%     of the domestic cat. Journal of the Association for Research in
%     Otolaryngology, 12(1):71--88, 2011. [1]http ]
%     
%     References
%     
%     1. http://link.springer.com/article/10.1007/s10162-010-0232-5
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2016_comparisonprocess.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA M-Signal M-Stats O-Statistics
%   #Author: Robert Baumgartner (2016), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



% Unbiased statistical index of rate discrimination acc. to Eq. 1 of
% may1997
sigma=zeros(size(tem.m,2),size(tar,2),size(tem.m,3),size(tem.m,4),size(tar,5)); % init
for itime = 1:size(tar.m,5)
  for itang = 1:size(tar.m,2)
    isd = repmat(tar.m(:,itang,:,:,itime),[1,size(tem.m,2),1,1]) - tem.m;
    sd = sqrt( repmat(tar.sd(:,itang,:,:,itime),[1,size(tem.sd,2),1,1]).^2 + tem.sd.^2);
    dprime = isd./sd;
    sigma(:,itang,:,:,itime) = nanmean(abs(dprime),1);
  end
end

end


