function sibin = baumgartner2014_binauralweighting(simon,varargin)
%BAUMGARTNER2014_BINAURALWEIGHTING Binaural combination of monaural similarity estimates
%   Usage:     sibin = baumgartner2014_binauralweighting(simon)
%
%   Input parameters:
%     simon   : monaural similarity indices
%
%   Output parameters:
%     sibin   : monaural similarity indices
%
%   BAUMGARTNER2014_BINAURALWEIGHTING(...) combines the monaural
%   similarity indices to binaural similartiy indices while accounting for
%   ipsilateral predominance.
%
%   BAUMGARTNER2014_BINAURALWEIGHTING accepts the following optional parameters:
%
%     'bwcoef',bwc   Set the binaural weighting coefficient bwc.
%                    Default value is 13 degrees.
%
%     'lat',lat      Set the apparent lateral angle of the target sound to
%                    lat. Default value is 0 degree (median SP).
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2014_binauralweighting.php


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

%% Binaural weighting, Eq.(6)

if size(simon,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.bwcoef)); % weight of left ear signal with 0 <= binw <= 1
    sibin = binw * simon(:,:,1,:,:) + (1-binw) * simon(:,:,2,:,:);
end

end


