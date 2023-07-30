function sigma = baumgartner2014_comparisonprocess(tar,tem)
%BAUMGARTNER2014_COMPARISONPROCESS Comparison with direction-specific templates
%   Usage:     sigma = baumgartner2014_comparisonprocess(tar,tem)
%
%   Input parameters:
%     tar     : internal representation(s) of target profile(s)
%     tem     : internal templates
%
%   Output parameters:
%     sigma   : internal distance metric
%
%   BAUMGARTNER2014_COMPARISONPROCESS(...) performs spectro-spatial
%   comparison process between internal representation of incoming sound
%   (target) and the templates of the sagittal plane
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2014_comparisonprocess.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA CircStat M-SIGNAL M-Stats O-Statistics
%   #Author: Robert Baumgartner (2014), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Comparison process, Eq.(4)
sigma=zeros(size(tem,2),size(tar,2),size(tem,3),size(tem,4),size(tar,5)); % init
for itime = 1:size(tar,5)
  for itang = 1:size(tar,2)
    isd = repmat(tar(:,itang,:,:,itime),[1,size(tem,2),1,1]) - tem;
    sigma(:,itang,:,:,itime) = mean(abs(isd),1);
  end
end

end


