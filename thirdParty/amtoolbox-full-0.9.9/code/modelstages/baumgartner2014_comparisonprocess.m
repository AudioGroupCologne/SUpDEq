function sigma = baumgartner2014_comparisonprocess(tar,tem)
%baumgartner2014_comparisonprocess - Comparison with direction-specific templates
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
%     Acoustical Society of America, 136(2):791-802, 2014.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2014_comparisonprocess.php

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
sigma=zeros(size(tem,2),size(tar,2),size(tem,3),size(tem,4),size(tar,5)); % init
for itime = 1:size(tar,5)
  for itang = 1:size(tar,2)
    isd = repmat(tar(:,itang,:,:,itime),[1,size(tem,2),1,1]) - tem;
    sigma(:,itang,:,:,itime) = mean(abs(isd),1);
  end
end

end
