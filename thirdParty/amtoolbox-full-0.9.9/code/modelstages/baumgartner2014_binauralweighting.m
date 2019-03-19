function sibin = baumgartner2014_binauralweighting(simon,varargin)
%baumgartner2014_binauralweighting - Binaural combination of monaural similarity estimates
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
%     Acoustical Society of America, 136(2):791-802, 2014.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2014_binauralweighting.php

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

%% Binaural weighting, Eq.(6)

if size(simon,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.bwcoef)); % weight of left ear signal with 0 <= binw <= 1
    sibin = binw * simon(:,:,1) + (1-binw) * simon(:,:,2);
end

end
