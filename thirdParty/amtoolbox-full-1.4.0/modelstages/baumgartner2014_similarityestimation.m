function si = baumgartner2014_similarityestimation(sigma,varargin)
%BAUMGARTNER2014_SIMILARITYESTIMATION Similarity estimation with listener-specific sensitivity
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
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2014_similarityestimation.php


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


