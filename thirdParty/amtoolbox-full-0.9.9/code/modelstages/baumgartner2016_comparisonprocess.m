function sigma = baumgartner2016_comparisonprocess(tar,tem)
%baumgartner2016_comparisonprocess - Comparison with direction-specific templates
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
%     J. Acoust. Soc. Am., 101:2705-2719, 1997.
%     
%     L. A. J. Reiss, R. Ramachandran, and B. J. May. Effects of signal level
%     and background noise on spectral representations in the auditory nerve
%     of the domestic cat. Journal of the Association for Research in
%     Otolaryngology, 12(1):71-88, 2011. [1]http ]
%     
%     References
%     
%     1. http://link.springer.com/article/10.1007/s10162-010-0232-5
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2016_comparisonprocess.php

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

% AUTHOR: Robert Baumgartner

%% Comparison process, Eq.(4)

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
