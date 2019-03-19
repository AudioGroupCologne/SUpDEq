function si = baumgartner2014_similarityestimation(sigma,varargin)
%baumgartner2014_similarityestimation - Similarity estimation with listener-specific sensitivity
%   Usage:     si = baumgartner2014_similarityestimation(sigma)
%
%   Input parameters:
%     sigma   : internal distance metrics
%
%   Output parameters:
%     si      : similarity indices
%
%   BAUMGARTNER2014_SIMILARITYESTIMATION(...) maps internal distance
%   metrics to similarity indices according to listener-specific
%   sensitivity
%
%   BAUMGARTNER2014_SIMILARITYESTIMATION accepts the following optional parameters:
%
%     'S',S          Set the listener-specific sensitivity threshold 
%                    (threshold of the sigmoid link function representing 
%                    the psychometric link between transformation from the
%                    distance metric and similarity index) to S. 
%                    Default value is 1.
%
%     'gamma',G      Set the degree of selectivity 
%                    (slope of the sigmoid link function representing 
%                    the psychometric link between transformation from the
%                    distance metric and similarity index) to G. 
%                    Default value is 6.
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791-802, 2014.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2014_similarityestimation.php

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

definput.import={'baumgartner2014'};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% Similarity estimation, Eq.(5)

si=zeros(size(sigma)); % init
for ch = 1:size(sigma,3)
  for it = 1:size(sigma,2)
    si(:,it,ch) = 1+eps - (1+exp(-kv.gamma*(sigma(:,it,ch)-kv.S))).^-1;
  end
end

end
